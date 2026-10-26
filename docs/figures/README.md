# Coupling figures

These six schematics use one visual language: blue for the 2-D surface, teal for the 1-D drainage network, and ochre for exchange. White backgrounds, flat fills and restrained linework keep attention on the physical relationships.

| Role | Colour |
|---|---|
| Surface water | `#306B9A` |
| Drainage network | `#247D74` |
| Exchange | `#BE7137` |
| Main text | `#26333B` |

Primary labels and equations use 14–22 pt type on a 10-inch-wide source canvas; brief secondary notes use 13–13.5 pt. Equations are typeset with MathText. SVG glyphs are exported as vector paths, so formulas and fonts do not depend on the reader's installed fonts. Each SVG has a matching 300-dpi PNG for presentations and raster previews. These are conceptual schematics, not geometrically scaled sections.

## Regenerate

The editable source is [tools/draw_figures.py](../../tools/draw_figures.py). Install Matplotlib and run from the repository root:

```bash
python tools/draw_figures.py
```

To render only one figure, use `--only manhole_exchange`; to avoid replacing the published figures, use `--output build/figure-preview`. No network service, external font download or LaTeX installation is needed. This style is inspired by clean scientific-journal graphics, not a claim of journal submission compliance.

## Surface and network

![Surface and drainage exchange](surface.svg)

## Architecture

![Model architecture](layers.svg)

## Exchange cycle

The three panels separate surface-step exchange, the network advance at an exchange time, and fixed-interval return scheduling. Direct manhole surcharge is part of the signed manhole exchange; outfall and node-flooding returns enter the return scheduler. Reservoir symbols show water pending transfer, and the footer includes both interface stores in the coupled budget. The [methodology](../coupling_methodology.md#3-one-coupling-cycle) defines the signs, source slots and units. Surface flux-kernel placement is solver-specific and is not drawn.

![Coupling cycle](exchange_cycle.svg)

## Manhole regimes

The free-inflow expression is algebraically identical to the implemented weir law. The three branches are tested in order, with the same head tolerance as the model. Full symbol definitions remain in the [methodology](../coupling_methodology.md#5-manhole-exchange).

![Manhole exchange regimes and supply limits](manhole_exchange.svg)

## Finite receiving cell

The donor constraint includes reserved water, consistent with the methodology. The drawn conduit is an interior-node-to-outfall link, whose end-storage weights are 1 and 0 respectively.

![Finite-volume outfall exchange](outfall_cv.svg)

## Time synchronisation

The three panels show fixed-interval return delay, every-step exchange, and CFL-bounded target alignment. Bar lengths represent simulated duration, not wall-clock execution time. SWMM routes an interval when the surface reaches its right boundary; the two solvers execute sequentially. Equal surface bars illustrate a constant CFL limit within a target segment; the model recomputes that limit each step. The orange rate bars show release from the pending queue over the following interval, while the arrows in step mode show current-step return application.

![Exchange timing](time_sync.svg)
