#pragma once

#include <cmath>

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

// Divide the remaining interval into equal steps no longer than the CFL limit.
inline double coupling_equal_substep(double dt, double remain)
{
	if (!(remain > 0.0)) return 0.0;
	if (!(dt > 0.0)) return dt;
	const double tol = 1.0e-9 * (remain > 1.0 ? remain : 1.0);
	if (dt >= remain - tol) return remain;
	double k = std::ceil(remain / dt - 1.0e-9);
	if (k < 1.0) k = 1.0;
	return remain / k;
}

// Land on the next exchange or end time using equal substeps.
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
	if (dt < 0.0) dt = 0.0;
	return coupling_equal_substep(dt, remain);
}

// Fixed-interval positive returns are spread over the next interval.
// Reverse outfall flow is immediate; the final step releases any remainder.
struct SwmmReturnSpread
{
	int count;          // coupling source slots (surface cells, then outfall cells)
	double* pending;    // m3 still to be applied
	double* rate;       // source per unit cell width (m2/s) during the current interval
};

inline void swmm_return_init(SwmmReturnSpread* s, int count)
{
	s->count = count > 0 ? count : 0;
	s->pending = s->count ? new double[s->count]() : 0;
	s->rate = s->count ? new double[s->count]() : 0;
}

inline void swmm_return_free(SwmmReturnSpread* s)
{
	delete[] s->pending; delete[] s->rate;
	s->pending = 0; s->rate = 0; s->count = 0;
}

inline double swmm_return_pending(const SwmmReturnSpread* s)
{
	double v = 0.0;
	for (int i = 0; i < s->count; ++i) v += s->pending[i];
	return v;
}

// Convert interval-mean sources (m2/s) into immediate withdrawals and queued returns.
template<typename T>
inline void swmm_return_schedule(SwmmReturnSpread* s, T* values,
	double exchange_seconds, double next_seconds, double width, double step_seconds)
{
	for (int i = 0; i < s->count; ++i)
	{
		const double q = (double)values[i];
		if (q < 0.0) values[i] = (T)(q * exchange_seconds / step_seconds);
		else { s->pending[i] += q * width * exchange_seconds; values[i] = (T)0; }
		if (next_seconds > 1.0e-9)
			s->rate[i] = s->pending[i] > 0.0 ? s->pending[i] / (next_seconds * width) : 0.0;
		else
		{
			if (s->pending[i] > 0.0) values[i] += (T)(s->pending[i] / (step_seconds * width));
			s->pending[i] = 0.0; s->rate[i] = 0.0;
		}
	}
}

// Every 2D step of an interval: adds the return source to values[] and returns
// the volume applied (m3) for the audit.
template<typename T>
inline double swmm_return_apply(SwmmReturnSpread* s, T* values, double width, double dt)
{
	double total = 0.0;
	if (!(dt > 0.0)) return 0.0;
	for (int i = 0; i < s->count; ++i)
	{
		if (!(s->rate[i] > 0.0) || !(s->pending[i] > 0.0)) continue;
		T r = (T)s->rate[i];
		double v = (double)r * width * dt;
		if (v > s->pending[i]) { r = (T)(s->pending[i] / (width * dt)); v = (double)r * width * dt; }
		values[i] += r;
		s->pending[i] -= v;
		total += v;
	}
	return total;
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
