#pragma once

// SWMM exchange clock. every_step keeps one exchange per 2D step.
// A positive interval accumulates 2D steps until their durations sum to
// that many seconds, then exchanges once.
struct SwmmCoupleClock
{
	int every_step;
	double interval;
	double elapsed;
	double* volume;
	int count;
};

inline void swmm_couple_clock_init
(
	SwmmCoupleClock* clock,
	int every_step,
	double interval,
	int count
)
{
	clock->every_step = every_step;
	clock->interval = interval;
	clock->elapsed = 0.0;
	clock->count = count;
	clock->volume = 0;
	if (!every_step && count > 0)
	{
		clock->volume = new double[count]();
	}
}

inline void swmm_couple_clock_free(SwmmCoupleClock* clock)
{
	delete[] clock->volume;
	clock->volume = 0;
}

// Shorten dt so the accumulated exchange time lands on the interval,
// the 2D end time, or the SWMM end time.
inline double swmm_couple_clip_dt
(
	const SwmmCoupleClock* clock,
	double t,
	double dt,
	double sim_end,
	double swmm_end
)
{
	if (clock->every_step) return dt;
	double remain = clock->interval - clock->elapsed;
	const double remain_sim = sim_end - t;
	const double remain_swmm = swmm_end - t;
	if (remain_sim < remain) remain = remain_sim;
	if (remain_swmm < remain) remain = remain_swmm;
	if (dt > remain) dt = remain;
	if (dt < 0.0) dt = 0.0;
	return dt;
}

inline void swmm_couple_accumulate
(
	SwmmCoupleClock* clock,
	Vconnect* verlink,
	double dt
)
{
	if (clock->every_step || dt <= 0.0) return;
	for (int i = 0; i < verlink->count_index; ++i)
	{
		clock->volume[i] += verlink->flow_1D_node[i] * dt;
	}
	clock->elapsed += dt;
}

inline int swmm_couple_due
(
	const SwmmCoupleClock* clock,
	double t_after,
	double sim_end,
	double swmm_end
)
{
	if (clock->every_step) return 1;
	if (clock->elapsed <= 0.0) return 0;
	const double eps = 1.0e-6;
	if (clock->elapsed + eps >= clock->interval) return 1;
	if (t_after + eps >= sim_end) return 1;
	if (t_after + eps >= swmm_end) return 1;
	return 0;
}

inline void swmm_couple_write_inflow
(
	SwmmCoupleClock* clock,
	Vconnect* verlink
)
{
	for (int i = 0; i < Nobjects[NODE]; ++i) Node_Inflow[i] = 0.0;
	if (clock->elapsed <= 0.0) return;
	for (int i = 0; i < verlink->count_index; ++i)
	{
		Node_Inflow[verlink->node_Vindex[i]] +=
				clock->volume[i] / clock->elapsed;
	}
}

inline void swmm_couple_advance(SwmmCoupleClock* clock)
{
	if (clock->every_step) return;
	for (int i = 0; i < clock->count; ++i) clock->volume[i] = 0.0;
	clock->elapsed = 0.0;
}
