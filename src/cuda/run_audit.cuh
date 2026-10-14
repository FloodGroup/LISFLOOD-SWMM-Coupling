#pragma once
#include <chrono>
#include <cstdio>
#include <string>
#include <stdexcept>
#include "../swmm5.h"
#include "cuda_stats.cuh"

// Lightweight interval diagnostics, not a change to the hydraulic equations.
class V11RunAudit {
    FILE* log_ = nullptr;
    FILE* water_ = nullptr;
    std::string status_;
    bool coupled_, finished_ = false;
    double initial2_, initial1_ = 0, expected_exchange_ = 0;
    std::chrono::steady_clock::time_point started_ = std::chrono::steady_clock::now();
    void status(const char* state, double t) {
        FILE* f=fopen(status_.c_str(),"w");
        if (!f) throw std::runtime_error("V11 cannot write run status: "+status_);
        fprintf(f,"engine=V11\nstate=%s\nsimulation_time_s=%.12g\nwall_since_solver_setup_s=%.3f\n",state,t,wall());
        int bad=ferror(f); if(fclose(f)!=0) bad=1;
        if(bad) throw std::runtime_error("V11 status write failed: "+status_);
    }
public:
    V11RunAudit(const char* root, const char* solver, bool coupled, double initial2, double target)
        : status_(std::string(root)+".status"),coupled_(coupled),initial2_(initial2) {
        log_=fopen((std::string(root)+".audit.log").c_str(),"w");
        water_=fopen((std::string(root)+".water_balance.csv").c_str(),"w");
        if(!log_ || !water_) throw std::runtime_error("V11 cannot create audit files");
        double s[5]={}; if(coupled_ && swmm_getCoupledAudit(s)==0) initial1_=s[0];
        fprintf(log_,"engine=V11 solver=%s precision=%s target_s=%.12g initial_2d_m3=%.12g initial_1d_m3=%.12g\n",solver,NUMERIC_TYPE_NAME,target,initial2_,initial1_);
        fprintf(log_,"Net-rain positive and negative volumes are NOT gross rainfall and gross infiltration. Other_2d includes external boundaries and uniform rain. Coupling entries are signed from 2D perspective. Pending exchange is included in combined balance.\n");
        fprintf(water_,"time_2d_s,time_swmm_s,wall_s,dt_s,steps_2d,steps_swmm,nonconverged_swmm,storage_2d_m3,storage_1d_m3,pending_to_1d_m3,netrain_added_m3,netrain_removed_m3,nodata_drained_m3,inlet_to_2d_m3,inlet_to_1d_m3,outfall_to_2d_m3,outfall_to_1d_m3,other_2d_net_m3,swmm_all_sources_net_m3,expected_coupling_to_2d_m3,actual_coupling_to_2d_m3,coupling_application_gap_m3,combined_balance_residual_m3\n");
        fflush(log_); status("RUNNING",0);
    }
    double wall() const { return std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count(); }
    void exchange(double volume_to_2d) { expected_exchange_+=volume_to_2d; }
    void sample(const lis::StatsEntry& e, lis::cuda::StatsCollector& stats, double pending) {
        double a[9]={},s[5]={}; stats.audit_totals(a);
        if(coupled_) swmm_getCoupledAudit(s);
        double coupling=a[5]-a[6]+a[7]-a[8];
        double net2=a[0]-a[1];
        double other=net2-coupling-(a[2]-a[3])+a[4];
        double residual=(double)e.volume-initial2_+s[0]-initial1_+pending-net2-s[1];
        fprintf(water_,"%.12g,%.12g,%.3f,%.9g,%ld,%.0f,%.0f,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            (double)e.t,s[2],wall(),(double)e.dt,e.iteration,s[3],s[4],(double)e.volume,s[0],pending,
            a[2],a[3],a[4],a[5],a[6],a[7],a[8],other,s[1],expected_exchange_,coupling,coupling-expected_exchange_,residual);
        fprintf(log_,"PROGRESS t_s=%.9g wall_s=%.3f dt_s=%.9g swmm_steps=%.0f nonconverged=%.0f combined_residual_m3=%.9g coupling_gap_m3=%.9g\n",
            (double)e.t,wall(),(double)e.dt,s[3],s[4],residual,coupling-expected_exchange_);
        if(fflush(water_)!=0 || fflush(log_)!=0) throw std::runtime_error("V11 audit write failed");
    }
    void finish(double t,double target,int error) {
        const char* state=error ? "FAILED" : (t+1.e-4>=target ? "COMPLETED" : "PARTIAL");
        fprintf(log_,"END state=%s t_s=%.12g target_s=%.12g error=%d wall_s=%.3f\n",state,t,target,error,wall());
        if(fflush(log_)!=0) throw std::runtime_error("V11 final audit write failed");
        status(state,t); finished_=true;
        printf("V11 END state=%s time_s=%.9g target_s=%.9g error=%d wall_since_solver_setup_s=%.3f\n",state,t,target,error,wall());
    }
    ~V11RunAudit() {
        if(log_) fclose(log_);
        if(water_) fclose(water_);
        // An interrupted process leaves RUNNING, never a false COMPLETED marker.
    }
};
