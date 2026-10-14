#pragma once
#include "../lisflood.h"
#include "../geometry.h"
#include "../lisflood2/DataTypes.h"
#include <cmath>
#include <algorithm>
#include <limits>

template<class Allocator = std::allocator<NUMERIC_TYPE>>
class DynamicRain
{
public:
    DynamicRain
    (
        const char* filename,
        int verbose,
        const Allocator& = Allocator()
    );

    /**
     * Check rainfall data has same origin as the DEM.
     */
    bool has_same_origin
    (
        Pars* Parptr
    );

    /**
     * Check rainfall tile dimensions are an integer multiple of
     * the DEM cell dimensions.
     */
    bool is_tile_size_multiple_of_grid
    (
        Pars* Parptr
    );

    /**
     * Update the current solver time, which loads new rainfall
     * data if necessary.
     */
    void update_time
    (
        NUMERIC_TYPE t
    );

/**
     * Update the current solver time, which loads new rainfall
     * data if necessary. For GGM returns the bool yes if updated
     */
    bool update_time_SGC
    (
        NUMERIC_TYPE t
    );

    /**
     * The rainfall rate (m/s) at cell (i, j).
     * If the cell is outside the rainfall domain, the rainfall rate is zero.
     */
    NUMERIC_TYPE rate_at_cell
    (
        Pars* Parptr,
        int i,
        int j
    );

    NUMERIC_TYPE rate_at_cell_SGC
    (
        int i,
        int j,
        const NUMERIC_TYPE dx_col,
        const NUMERIC_TYPE dy_col,
        NUMERIC_TYPE tly
    );

    void update_H
    (
        Pars* Parptr,
        Solver *Solverptr,
        Arrays *Arrptr
    );

    void update_rain_grid_SGM
    (
        const NUMERIC_TYPE curr_time,
        NUMERIC_TYPE* rain_grid,
        const NUMERIC_TYPE* dem_grid,
        WetDryRowBound* wet_dry_bounds,
        const NUMERIC_TYPE* dx_col,
        const NUMERIC_TYPE* dy_col,
        const NUMERIC_TYPE tyl,
        const int grid_rows,
        const int grid_cols_padded,
        const NUMERIC_TYPE* cell_area_col
    );

    inline bool enabled()
    {
        return enabled_;
    }

    // Hold each interval's rate constant, but never integrate it past its end.
    NUMERIC_TYPE clip_step(NUMERIC_TYPE t, NUMERIC_TYPE step) const
    {
#if _NETCDF == 1
        if (enabled_) {
            for (size_t i=0;i<netcdf_.tlen;++i) {
                // Round the timestamp to the solver's representable seconds,
                // exactly as the frame selector does (avoid zero-progress steps).
                NUMERIC_TYPE boundary=netcdf_.times[i]*C(3600.0);
                if (boundary>(double)t) return (NUMERIC_TYPE)(std::min)((double)step,boundary-(double)t);
            }
        }
#endif
        return step;
    }

    void validate_end(NUMERIC_TYPE end) const
    {
#if _NETCDF == 1
        if (enabled_) {
            const size_t n=netcdf_.tlen;
            double coverage=(2.0*(double)netcdf_.times[n-1]-netcdf_.times[n-2])*3600.0;
            if ((double)end>coverage+1.e-4) {
                fprintf(stderr,"V11 rainfall coverage ends at %.9g s, requested %.9g s: %s\n",coverage,(double)end,netcdf_.source_name);
                exit(EXIT_FAILURE);
            }
        }
#endif
    }

	inline const Geometry& geometry()
    {
        return geometry_;
    }

    /**
     * Spatial map of rainfall rate (m/s) for the current time.
     */
    inline NUMERIC_TYPE* data()
    {
        return netcdf_.data;
    }

    ~DynamicRain();

private:
    bool enabled_ = false;
    NetCDFVariable netcdf_ = {};
    Geometry geometry_ = {};
    Allocator allocator_;
};
