#ifndef COUPLED_BOUNDARY_H
#define COUPLED_BOUNDARY_H
/* Internal units: feet, cubic feet, seconds. Public interface uses SI. */
int coupled_open(void);
void coupled_close(void);
int coupled_active(void);
int coupled_node(int node);
int coupled_outfall(int node);
int coupled_set_cell(int node, int cell, double bed, double depth, double area);
void coupled_begin_step(void);
void coupled_limit_flows(double dt);
double coupled_solve_node(int node, double dt);
double coupled_solve_cells(double dt);
int coupled_commit(double dt);
void coupled_account(double dt);
void coupled_begin_interval(void);
double coupled_exchange(int node);
double coupled_standard_exchange(int node);
void coupled_sync_inlets(void);
void coupled_si_inlets(int enabled);
double coupled_inlet_factor(void);
int coupled_ready(void);
double coupled_available(int node);
double coupled_old_volume(int node);
void coupled_remaining_interval(double seconds);
int coupled_restart(void* stream, int reading);
#endif
