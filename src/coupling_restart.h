#pragma once
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "lisflood.h"
#include "geometry.h"
#include "swmm5.h"
#include "swmm_couple.h"

// A restart is tied to the same solver, backend, numeric precision and mesh.
// Host fields include all prognostic variables, not just output snapshots.
struct CouplingRestartField {
    const char* name;
    NUMERIC_TYPE* data;
    std::size_t count;
};

inline void restart_bytes(FILE* f, void* p, std::size_t bytes, bool read) {
    if (bytes && (read ? std::fread(p, 1, bytes, f) : std::fwrite(p, 1, bytes, f)) != bytes)
        throw std::runtime_error("Incomplete coupled restart file");
}

template<class T> inline void restart_value(FILE* f, T& v, bool read) {
    restart_bytes(f, &v, sizeof(v), read);
}

inline std::uint64_t coupling_restart_mapping(const Vconnect* v) {
    std::uint64_t hash = 1469598103934665603ULL;
    auto add = [&hash](int value) { hash ^= (std::uint32_t)value; hash *= 1099511628211ULL; };
    add(v->count_index); add(v->count_dem_index); add(v->outfall_count);
    for (int i=0; i<v->count_index; ++i) { add(v->node_Vindex[i]); add(v->dem_Vindex[i]); }
    for (int i=0; i<v->outfall_count; ++i) { add(v->outfall_node_Vindex[i]); add(v->outfall_dem_Vindex[i]); }
    return hash;
}

inline double coupling_pending_in(const SwmmCoupleClock& clock) {
    double sum = 0.0;
    if (clock.volume) for (int i=0; i<clock.count; ++i) sum += clock.volume[i];
    return sum;
}

inline double coupling_step(const SwmmCoupleClock& clock, bool coupled,
    double time, double cfl_dt, double end, double swmm_now, double swmm_end,
    double rain_remain, double checkpoint_time)
{
    double remain=end-time;
    if(coupled) {
        remain=(std::min)(remain,swmm_end-(clock.every_step?swmm_now:time));
        if(!clock.every_step) remain=(std::min)(remain,clock.interval-clock.elapsed);
    }
    if(rain_remain<remain-RAIN_FRAME_TIME_TOL) remain=rain_remain;
    if(checkpoint_time>time+1.e-7) remain=(std::min)(remain,checkpoint_time-time);
    return coupling_equal_substep(cfl_dt,remain);
}

inline std::vector<CouplingRestartField> cpu_restart_fields(Pars& pars, Arrays& arrays, bool fv1)
{
    const std::size_t cells=(std::size_t)pars.xsz*pars.ysz;
    const std::size_t faces=(std::size_t)(pars.xsz+1)*(pars.ysz+1);
    std::vector<CouplingRestartField> fields={{"H",arrays.H,cells},{"maxH",arrays.maxH,cells},
        {"maxHtm",arrays.maxHtm,cells},{"initHtm",arrays.initHtm,cells},{"totalHtm",arrays.totalHtm,cells}};
    if(fv1) { fields.push_back({"HU",arrays.HU,cells}); fields.push_back({"HV",arrays.HV,cells}); }
    else { fields.push_back({"Qx",arrays.Qx,faces}); fields.push_back({"Qy",arrays.Qy,faces});
        fields.push_back({"Qxold",arrays.Qxold,faces}); fields.push_back({"Qyold",arrays.Qyold,faces}); }
    return fields;
}

// Full hydraulic restart. SWMM's native warm start remains a separate legacy
// option; this format also retains the interface reservoirs and time history.
inline void coupling_restart_io(const char* path, bool read, int backend,
    Pars& pars, Solver& solver, const Vconnect* mapping, bool coupled,
    SwmmCoupleClock& clock, SwmmReturnSpread& returns, double& exact_time,
    const std::vector<CouplingRestartField>& fields)
{
    if (coupled && !swmm_hydraulicRestartSupported())
        throw std::runtime_error("Coupled hydraulic restart requires no SWMM runoff, RDII, storage exfiltration, routing-interface input, events or active water-quality routing; continuous runs remain supported");
    FILE* f = std::fopen(path, read ? "rb" : "wb");
    if (!f) throw std::runtime_error(std::string("Cannot open coupled restart: ") + path);
    try {
        char magic[8] = {'L','F','S','W','M','H','0','1'};
        restart_bytes(f, magic, 8, read);
        if (std::memcmp(magic, "LFSWMH01", 8)) throw std::runtime_error("Invalid coupled restart signature");
        int layout[8] = {1, (int)sizeof(NUMERIC_TYPE), backend, pars.xsz, pars.ysz,
                         coupled ? 1 : 0, clock.count, returns.count};
        int stored[8]; std::memcpy(stored, layout, sizeof(layout));
        restart_bytes(f, stored, sizeof(stored), read);
        if (std::memcmp(stored, layout, sizeof(layout)))
            throw std::runtime_error("Restart solver/backend/precision/grid/source counts do not match");
        double geometry[4] = {(double)pars.blx, (double)pars.bly, (double)pars.dx, (double)pars.dy};
        double saved_geometry[4]; std::memcpy(saved_geometry, geometry, sizeof(geometry));
        restart_bytes(f, saved_geometry, sizeof(saved_geometry), read);
        if (std::memcmp(saved_geometry, geometry, sizeof(geometry))) throw std::runtime_error("Restart grid geometry does not match");
        std::uint64_t signature = coupled ? coupling_restart_mapping(mapping) : 0, saved_signature = signature;
        restart_value(f, saved_signature, read);
        if (saved_signature != signature) throw std::runtime_error("Restart cell/node mapping does not match");
        int mode=clock.every_step; double interval=clock.interval;
        restart_value(f, mode, read); restart_value(f, interval, read);
        if (mode!=clock.every_step || (!mode && std::fabs(interval-clock.interval)>1.e-12))
            throw std::runtime_error("Restart exchange interval does not match");
        restart_value(f, exact_time, read);
        if (!std::isfinite(exact_time) || exact_time < 0 || exact_time >= (double)solver.Sim_Time)
            if (read) throw std::runtime_error("Restart time must precede sim_time; use the original absolute forcing timeline");
        restart_value(f, clock.elapsed, read);
        restart_value(f, solver.Tstep, read); restart_value(f, solver.MinTstep, read);
        restart_value(f, solver.itCount, read);
        restart_value(f, pars.SaveNo, read); restart_value(f, pars.SaveTotal, read); restart_value(f, pars.MassTotal, read);
        if (clock.volume) restart_bytes(f, clock.volume, clock.count*sizeof(double), read);
        restart_bytes(f, returns.pending, returns.count*sizeof(double), read);
        restart_bytes(f, returns.rate, returns.count*sizeof(double), read);
        for (int i=0; i<returns.count; ++i)
            if (!std::isfinite(returns.pending[i]) || returns.pending[i]<-1.e-8 ||
                !std::isfinite(returns.rate[i]) || returns.rate[i]<0)
                throw std::runtime_error("Invalid pending return reservoir in restart");
        std::uint64_t count=fields.size(); restart_value(f, count, read);
        if (count!=fields.size()) throw std::runtime_error("Restart surface-field layout does not match");
        for (const auto& field: fields) {
            char name[40] = {}; std::strncpy(name, field.name, sizeof(name)-1);
            char saved_name[40]; std::memcpy(saved_name,name,sizeof(name));
            std::uint64_t size=field.count;
            restart_bytes(f,saved_name,sizeof(saved_name),read); restart_value(f,size,read);
            if (std::memcmp(name,saved_name,sizeof(name)) || size!=field.count || !field.data)
                throw std::runtime_error("Restart surface-field name or size does not match");
            restart_bytes(f,field.data,size*sizeof(NUMERIC_TYPE),read);
            if (read) for (std::size_t i=0;i<size;++i)
                if (!std::isfinite((double)field.data[i])) throw std::runtime_error("Nonfinite restart surface field");
        }
        if (coupled && swmm_hydraulicRestart(f,read?1:0))
            throw std::runtime_error("SWMM hydraulic restart is incompatible or incomplete");
        std::uint64_t footer=0x51ABCF2786D01394ULL, saved_footer=footer;
        restart_value(f,saved_footer,read);
        if (saved_footer!=footer || (read && std::fgetc(f)!=EOF)) throw std::runtime_error("Invalid restart footer");
        if (!read && std::fflush(f)) throw std::runtime_error("Restart flush failed");
        if (std::fclose(f)) { f=nullptr; throw std::runtime_error("Restart close failed"); }
        f=nullptr;
        solver.t=(NUMERIC_TYPE)exact_time;
        std::printf("RESTART %s time_s=%.12g pending_to_1d_m3=%.12g pending_to_2d_m3=%.12g exchange_elapsed_s=%.12g\n",
            read?"restored":"saved",exact_time,coupling_pending_in(clock),swmm_return_pending(&returns),clock.elapsed);
    } catch (...) { if(f) std::fclose(f); throw; }
}
