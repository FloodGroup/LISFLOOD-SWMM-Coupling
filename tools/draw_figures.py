"""Generate the documentation schematics as vector SVG and high-resolution PNG.

Requires matplotlib. Run from any directory with Python 3.
Colours encode model domains consistently; equations use MathText, not Unicode
approximations. SVG glyphs are paths so browsers need no locally installed fonts.
"""
from pathlib import Path
import argparse
from io import StringIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, Circle


BLUE = '#306B9A'
TEAL = '#247D74'
ORANGE = '#BE7137'
INK = '#26333B'
GRAY = '#596870'
LINE = '#CCD5DA'
PALE_BLUE = '#E6EFF6'
PALE_TEAL = '#E4F0ED'
PALE_ORANGE = '#F8EDE3'
SOIL = '#F0EFEB'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 15,
    'mathtext.fontset': 'dejavusans', 'svg.fonttype': 'path',
    'svg.hashsalt': 'lisflood-swmm-documentation',
    'savefig.facecolor': 'white', 'axes.unicode_minus': True,
})


def canvas(height, title, subtitle=None):
    fig = plt.figure(figsize=(10, height / 100), dpi=160)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, 1000), ylim=(height, 0))
    ax.axis('off')
    text(ax, 35, 31, title, size=21, weight='bold')
    if subtitle:
        text(ax, 35, 72, subtitle, size=13.5, color=GRAY)
    return fig, ax


def text(ax, x, y, s, size=15, color=INK, weight='normal', ha='left', va='center'):
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
                   ha=ha, va=va, linespacing=1.45, zorder=6)


def rect(ax, x, y, w, h, fill='white', edge=LINE, lw=1.1, **kwargs):
    if kwargs.pop('dashed', False):
        kwargs['linestyle'] = (0, (4, 3))
    p = Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge, linewidth=lw, **kwargs)
    ax.add_patch(p)
    return p


def line(ax, xs, ys, color=INK, lw=1.5, dashed=False):
    ax.plot(xs, ys, color=color, lw=lw, linestyle=(0, (4, 3)) if dashed else '-',
            solid_capstyle='round')


def arrow(ax, x1, y1, x2, y2, color=ORANGE, lw=2, both=False):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                 arrowstyle='<->' if both else '-|>', mutation_scale=16,
                 color=color, linewidth=lw, shrinkA=0, shrinkB=0, zorder=5))


def routed(ax, points, color=GRAY):
    line(ax, [p[0] for p in points[:-1]], [p[1] for p in points[:-1]], color)
    arrow(ax, *points[-2], *points[-1], color=color, lw=1.5)


def block(ax, x, y, w, h, title, body, color, pale):
    rect(ax, x, y, w, h, pale, edge='none')
    rect(ax, x, y, 4, h, color, edge='none')
    text(ax, x + 20, y + 28, title, size=15.5, color=color, weight='bold')
    text(ax, x + 20, y + (59 if '\n' in body else 66), body,
         size=13.5, va='top' if '\n' in body else 'center')


def surface():
    fig, ax = canvas(530, 'Two-way exchange between surface and drainage')
    text(ax, 40, 90, '2-D surface', color=BLUE, size=18, weight='bold')
    text(ax, 40, 130, 'LISFLOOD-FP', color=GRAY, size=14)
    ax.add_patch(Polygon([(40, 220), (660, 220), (735, 410), (960, 410),
                          (960, 455), (40, 455)], facecolor=SOIL, edgecolor='none'))
    rect(ax, 40, 195, 620, 25, PALE_BLUE, edge='none')
    line(ax, [40, 660], [195, 195], BLUE, 2)
    line(ax, [40, 660, 735, 960], [220, 220, 410, 410], INK, 2)
    rect(ax, 735, 195, 225, 215, PALE_BLUE, edge='none')
    line(ax, [735, 960], [195, 195], BLUE, 2)
    rect(ax, 735, 160, 225, 250, 'none', BLUE, dashed=True)
    rect(ax, 120, 350, 620, 32, PALE_TEAL, INK, 1.6, zorder=3)
    for x in (225, 485):
        rect(ax, x-18, 220, 36, 145, 'white', INK, 1.6, zorder=3)
        water_top = 295 if x == 225 else 222
        rect(ax, x-16, water_top, 32, 364-water_top, PALE_TEAL, edge='none', zorder=4)
    arrow(ax, 225, 200, 225, 322)
    arrow(ax, 485, 322, 485, 180)
    text(ax, 225, 165, 'Inlet drainage', size=15, ha='center', color=ORANGE)
    text(ax, 485, 165, 'Surcharge', size=15, ha='center', color=ORANGE)
    text(ax, 790, 95, 'Receiving cell', color=BLUE, weight='bold', size=16)
    text(ax, 790, 135, 'river or coast', color=GRAY, size=14)
    arrow(ax, 695, 336, 795, 336, both=True)
    text(ax, 790, 270, 'COUPLED', color=ORANGE, ha='center', weight='bold')
    text(ax, 790, 305, 'outfall', color=ORANGE, ha='center', size=14)
    text(ax, 350, 418, '1-D drainage network · SWMM', color=TEAL, ha='center', weight='bold')
    line(ax, [40, 960], [480, 480], LINE, 1)
    text(ax, 40, 506, 'Water levels drive exchange; available storage limits withdrawal.', size=14)
    return fig


def layers():
    fig, ax = canvas(550, 'Three components, one coupled calculation',
                     'CPU and CUDA surface solvers use the same exchange laws.')
    block(ax, 35, 140, 280, 180, '2-D surface',
          'LISFLOOD-FP\nACC or FV1\nCPU / CUDA', BLUE, PALE_BLUE)
    block(ax, 360, 140, 280, 180, 'Coupling',
          'Cell–node mapping\nExchange clock\nFlow and volume limits', ORANGE, PALE_ORANGE)
    block(ax, 685, 140, 280, 180, '1-D network',
          'SWMM dynamic wave\nAdaptive substeps\nCoupled continuity', TEAL, PALE_TEAL)
    arrow(ax, 320, 230, 355, 230, both=True, color=GRAY)
    arrow(ax, 645, 230, 680, 230, both=True, color=GRAY)
    text(ax, 500, 367, 'Surface heads + receiving-cell geometry + inlet volume',
         size=14, color=BLUE, ha='center')
    arrow(ax, 100, 395, 900, 395, color=BLUE)
    text(ax, 500, 436, 'Node heads + storage + outfall / flooding volumes',
         size=14, color=TEAL, ha='center')
    arrow(ax, 900, 463, 100, 463, color=TEAL)
    text(ax, 500, 518, 'CPU: direct depth access     |     CUDA: gather coupled-cell depths only',
         size=13.5, ha='center', color=GRAY)
    return fig


def _panel_heading(ax, letter, title, y, note=None):
    text(ax, 35, y, letter, size=18, weight='bold')
    text(ax, 68, y, title, size=16, weight='bold')
    if note:
        text(ax, 965, y, note, size=13, color=GRAY, ha='right')


def _flow_card(ax, x, y, w, h, title, color, pale, title_size=15):
    boundary = rect(ax, x, y, w, h, pale, edge=LINE, lw=0.8)
    boundary._constrain_text = True
    rect(ax, x, y, 4, h, color, edge='none')
    text(ax, x+18, y+29, title, size=title_size, color=color, weight='bold')


def _reservoir_symbol(ax, x, y, w, h, color, pale):
    rect(ax, x, y+h*0.38, w, h*0.62, pale, edge='none')
    line(ax, [x, x, x+w, x+w], [y, y+h, y+h, y], color, 1.4)
    line(ax, [x, x+w], [y+h*0.38, y+h*0.38], color, 1.4)


def _time_strip(ax, x, y, widths, color, pale, height=34):
    for w in widths:
        rect(ax, x, y, w, height, pale, edge=color, lw=0.95)
        x += w


def exchange_cycle():
    fig, ax = canvas(1180, 'Two-way coupling and water in transit',
                     'Surface forcing, interval routing and return scheduling are shown separately.')

    _panel_heading(ax, 'a', 'Every surface step', 120,
                   'Repeat until an exchange is due')
    _flow_card(ax, 35, 154, 280, 130, '1  Apply net forcing', BLUE, PALE_BLUE)
    text(ax, 53, 227, 'Apply signed net forcing', size=13.5)
    text(ax, 53, 259, 'before inlet withdrawal', size=13.5)
    _flow_card(ax, 360, 154, 280, 130, '2  Manhole exchange', ORANGE, PALE_ORANGE)
    text(ax, 378, 227, 'Use updated water levels', size=13.5)
    text(ax, 378, 259, 'Apply donor-volume limits', size=13.5)
    _flow_card(ax, 685, 154, 280, 130, '3  Sources + volume', BLUE, PALE_BLUE)
    text(ax, 703, 222, 'Manholes + return queue', size=13.5)
    text(ax, 825, 261, r'$V_i\leftarrow V_i+Q_i\,\Delta t$', size=16, ha='center')
    arrow(ax, 315, 219, 355, 219, color=GRAY, lw=1.6)
    arrow(ax, 640, 219, 680, 219, color=GRAY, lw=1.6)
    text(ax, 35, 316, r'$Q_i>0$: surface to network; $Q_i<0$: direct manhole surcharge.',
         size=13.5, color=GRAY)
    line(ax, [35, 965], [341, 341], LINE, 0.9)

    _panel_heading(ax, 'b', 'At an exchange time', 373,
                   r'$t_n=t_{n-1}+T_n$')
    text(ax, 35, 408, 'Supply the current receiving-cell bed, depth and area before routing.',
         size=13.5, color=GRAY)
    _flow_card(ax, 35, 435, 260, 149, 'Inlet volume to rate', ORANGE, PALE_ORANGE)
    text(ax, 165, 508, r'$\overline{Q}_i=\frac{V_i}{T_n}$', size=18, ha='center')
    text(ax, 165, 551, r'$T_n=t_n-t_{n-1}$', size=15, ha='center', color=GRAY)
    _flow_card(ax, 365, 435, 280, 149, 'Advance SWMM', TEAL, PALE_TEAL)
    text(ax, 505, 499, 'Adaptive routing steps', size=13.5, ha='center')
    _time_strip(ax, 388, 521, [45, 27, 63, 38, 61], TEAL, PALE_TEAL, height=19)
    text(ax, 505, 565, r'$\sum_k\delta t_k=T_n$', size=15, ha='center', color=TEAL)
    _flow_card(ax, 710, 435, 255, 149, 'Collect net returns', TEAL, PALE_TEAL)
    text(ax, 728, 503, 'Outfalls + node flooding', size=13.5)
    text(ax, 837, 550, r'$\Delta V_s$', size=22, color=ORANGE, ha='center')
    arrow(ax, 295, 510, 360, 510, color=ORANGE, lw=1.7)
    arrow(ax, 645, 510, 705, 510, color=TEAL, lw=1.7)
    line(ax, [35, 965], [608, 608], LINE, 0.9)

    _panel_heading(ax, 'c', 'Return scheduling in fixed-interval mode', 640,
                   'Per source slot')
    text(ax, 35, 675, r'$\Delta V_s$ is positive towards the surface; aggregate by source slot before scheduling.',
         size=13.2, color=GRAY)
    rect(ax, 447, 700, 106, 45, PALE_ORANGE, edge=ORANGE, lw=1.1)
    text(ax, 500, 723, r'$\Delta V_s$', size=19, color=ORANGE, ha='center')
    line(ax, [500, 500, 258], [745, 762, 762], ORANGE, 1.5)
    line(ax, [500, 742], [762, 762], ORANGE, 1.5)
    arrow(ax, 258, 762, 258, 790, color=ORANGE, lw=1.5)
    arrow(ax, 742, 762, 742, 790, color=ORANGE, lw=1.5)

    _flow_card(ax, 35, 795, 445, 185, 'Reverse outfall flow', BLUE, PALE_BLUE)
    text(ax, 462, 824, r'$\Delta V_s<0$', size=15, color=ORANGE, ha='right')
    _reservoir_symbol(ax, 68, 853, 62, 55, BLUE, 'white')
    arrow(ax, 99, 890, 99, 942, color=ORANGE, lw=1.7)
    text(ax, 99, 960, '2-D cell', size=13, color=BLUE, ha='center')
    text(ax, 304, 868, r'$\Delta h_s=\frac{\Delta V_s}{A_s}$', size=18, color=BLUE, ha='center')
    text(ax, 176, 916, 'Deduct now; no delayed sink.', size=13.2)
    text(ax, 176, 953, 'Use the available 2-D supply.', size=13.2)

    _flow_card(ax, 520, 795, 445, 185, 'Pending return reservoir', ORANGE, PALE_ORANGE)
    text(ax, 947, 824, r'$\Delta V_s\geq0$', size=15, color=ORANGE, ha='right')
    _reservoir_symbol(ax, 553, 853, 62, 55, ORANGE, 'white')
    text(ax, 584, 883, r'$P_s$', size=16, color=ORANGE, ha='center')
    arrow(ax, 584, 910, 584, 942, color=BLUE, lw=1.7)
    text(ax, 584, 960, 'to 2-D', size=13, color=BLUE, ha='center')
    text(ax, 791, 865, r'$P_s\leftarrow P_s+\Delta V_s$', size=17, ha='center')
    text(ax, 791, 911, r'$Q_s^{\mathrm{ret}}=\frac{P_s}{T_{\mathrm{next}}}$', size=17,
         color=ORANGE, ha='center')
    text(ax, 791, 953, 'Release over the next interval.', size=13.2, ha='center')

    text(ax, 35, 1016, 'Step mode: current-step return. Final exchange: release all remaining returns.',
         size=13.5, color=GRAY)
    line(ax, [35, 965], [1043, 1043], LINE, 0.9)
    text(ax, 500, 1073, 'Interface reservoirs belong to the coupled water budget',
         size=14, weight='bold', ha='center')
    text(ax, 500, 1113, r'$S_{\mathrm{total}}=S_{\mathrm{2D}}+S_{\mathrm{1D}}+P_{\to\mathrm{1D}}+P_{\to\mathrm{2D}}$',
         size=18, ha='center')
    text(ax, 500, 1152, 'Restart retains both reservoirs, return rates and exchange progress, with both model states.',
         size=13, color=GRAY, ha='center')
    return fig


def manhole_exchange():
    fig, ax = canvas(870, 'Manhole exchange: three hydraulic regimes',
                     r'$Q_*$ is the hydraulic demand, positive in the direction of each arrow.')
    panels = [
        (25, 'a  Surcharge', 'sewer → surface', 300, 245,
         r'$Q_*=C_0 A_0\sqrt{2g\,\Delta H}$', r'$\Delta H=H_{\mathrm{1D}}-H_{\mathrm{2D}}$',
         r'$H_{\mathrm{1D}}>Z+\varepsilon$', r'$H_{\mathrm{1D}}>H_{\mathrm{2D}}+\varepsilon$'),
        (350, 'b  Free inflow', 'surface → sewer', 245, 400,
         r'$Q_*=C_w w\sqrt{2g}\,h_c^{3/2}$', r'$h_c=H_{\mathrm{2D}}-Z$',
         r'$H_{\mathrm{1D}}\leq Z+\varepsilon$', r'$H_{\mathrm{2D}}>Z+\varepsilon$'),
        (675, 'c  Submerged inflow', 'surface → sewer', 245, 285,
         r'$Q_*=C_1 A_1\sqrt{2g\,\Delta H}$', r'$\Delta H=H_{\mathrm{2D}}-H_{\mathrm{1D}}$',
         r'$H_{\mathrm{1D}}\geq Z-\varepsilon$', r'$H_{\mathrm{2D}}>H_{\mathrm{1D}}+\varepsilon$'),
    ]
    for x, title, direction, h2, h1, eq, head, c1, c2 in panels:
        text(ax, x+8, 130, title, size=15.5, weight='bold')
        text(ax, x+8, 163, direction, size=14, color=ORANGE)
        rect(ax, x+8, 330, 284, 146, SOIL, edge='none')
        rect(ax, x+8, h2, 284, 330-h2, PALE_BLUE, edge='none')
        line(ax, [x+8, x+292], [h2, h2], BLUE, 2)
        rect(ax, x+25, 438, 250, 28, PALE_TEAL, INK)
        rect(ax, x+126, 330, 48, 122, 'white', INK)
        rect(ax, x+128, max(h1, 330), 44, 452-max(h1, 330), PALE_TEAL, edge='none')
        line(ax, [x+8, x+126], [330, 330], INK, 2)
        line(ax, [x+174, x+292], [330, 330], INK, 2)
        line(ax, [x+80, x+292], [h1, h1], TEAL, 1.5, dashed=True)
        text(ax, x+285, h1-15, r'$H_{\mathrm{1D}}$', size=15, color=TEAL, ha='right')
        text(ax, x+12, h2-15, r'$H_{\mathrm{2D}}$', size=15, color=BLUE)
        text(ax, x+12, 352, '$Z$', size=15)
        if title.startswith('a'):
            arrow(ax, x+150, 420, x+150, 275)
        else:
            arrow(ax, x+150, 270, x+150, 420)
        text(ax, x+150, 515, eq, size=15, ha='center')
        text(ax, x+150, 555, head, size=14, ha='center', color=GRAY)
        text(ax, x+150, 610, c1, size=14, ha='center')
        text(ax, x+150, 645, c2, size=14, ha='center')
    line(ax, [30, 970], [687, 687], LINE, 1)
    text(ax, 35, 720, 'Surcharge limit', size=15, color=TEAL, weight='bold')
    text(ax, 525, 720, 'Inlet limit', size=15, color=BLUE, weight='bold')
    text(ax, 35, 760, r'$Q_{\mathrm{out}}=\min\!\left(Q_*,Q_{\mathrm{out,max}},\frac{V_{\mathrm{avail}}}{\Delta t}\right)$', size=14)
    text(ax, 525, 760, r'$Q_{\mathrm{in}}=\min\!\left(Q_*,Q_{\mathrm{in,max}},\frac{h\,\Delta x\,\Delta y}{n\,\Delta t}\right)$', size=14)
    text(ax, 35, 802, r'$V_{\mathrm{avail}}=V_{\mathrm{stored}}+V_{\mathrm{pending}}$', size=13, color=GRAY)
    text(ax, 525, 802, '$n$: manholes sharing the cell', size=13, color=GRAY)
    text(ax, 500, 845, r'Branches evaluated a → b → c; $\varepsilon=10^{-3}$ m. Otherwise $Q=0$.',
         size=13.5, ha='center', color=GRAY)
    return fig


def outfall_cv():
    fig, ax = canvas(710, 'A COUPLED outfall connects to a finite receiving cell',
                     'The outfall is an exchange face, not a storage node.')
    text(ax, 65, 133, 'Interior node', color=TEAL, weight='bold')
    text(ax, 433, 133, 'Outfall face', color=ORANGE, weight='bold')
    text(ax, 715, 133, '2-D cell', color=BLUE, weight='bold')
    ax.add_patch(Polygon([(40, 200), (450, 200), (565, 405), (960, 405),
                          (960, 430), (40, 430)], facecolor=SOIL, edgecolor='none'))
    rect(ax, 565, 225, 395, 180, PALE_BLUE, edge='none')
    rect(ax, 565, 175, 395, 230, 'none', BLUE, dashed=True)
    line(ax, [565, 960], [225, 225], BLUE, 2)
    rect(ax, 110, 320, 460, 35, PALE_TEAL, INK, 1.5)
    rect(ax, 100, 195, 55, 160, 'white', TEAL, 1.5)
    rect(ax, 102, 265, 51, 88, PALE_TEAL, edge='none')
    line(ax, [102, 155], [265, 265], TEAL, 2)
    text(ax, 170, 265, r'$H_{\mathrm{1D}}$', color=TEAL, size=16)
    line(ax, [565, 565], [280, 385], ORANGE, 1.6, dashed=True)
    ax.add_patch(Circle((565, 337), 5, facecolor=ORANGE, edgecolor='none'))
    arrow(ax, 490, 305, 655, 305, both=True)
    text(ax, 650, 275, 'Two-way exchange', color=ORANGE, size=14)
    text(ax, 725, 204, r'$H_{\mathrm{ext}}$', color=BLUE, size=16)
    arrow(ax, 910, 230, 910, 400, color=BLUE, lw=1.4, both=True)
    text(ax, 887, 365, '$h$', color=BLUE, size=17)
    text(ax, 757, 393, r'$z_{\mathrm{bed}}$', color=INK, size=14)
    text(ax, 55, 460, 'Conduit storage → interior node', color=TEAL, size=14)
    text(ax, 565, 460, 'Zero storage weight at the face', color=ORANGE, size=14)
    line(ax, [35, 965], [490, 490], LINE, 1)
    text(ax, 45, 524, 'Receiving-cell stage', color=BLUE, weight='bold')
    text(ax, 535, 524, 'Integrated exchange', color=ORANGE, weight='bold')
    text(ax, 45, 569, r'$V_{\mathrm{cell}}=Ah,\quad H_{\mathrm{ext}}=z_{\mathrm{bed}}+\frac{V_{\mathrm{cell}}}{A}$', size=16)
    text(ax, 535, 569, r'$\Delta V_{\mathrm{ex}}=\sum_k Q_{\mathrm{net},k}\,\delta t_k$', size=17)
    text(ax, 45, 635, 'Finite donor supply', size=14, weight='bold')
    text(ax, 320, 635, r'$\sum Q_{\mathrm{out}}\leq\frac{V^n-V_{\mathrm{reserved}}}{\delta t}+\sum Q_{\mathrm{in}}$', size=17)
    text(ax, 45, 680, 'Reserved volume is kept for exchanges already promised to the surface.', size=13.5, color=GRAY)
    return fig


def tiles(ax, x, y, widths, color, pale, labels=False):
    for i, w in enumerate(widths):
        rect(ax, x, y, w, 42, pale, color, 1)
        if labels:
            text(ax, x+w/2, y+21, rf'$\Delta t_{{{i+1}}}$', size=14, ha='center', color=color)
        x += w


def time_sync():
    fig, ax = canvas(1090, 'Time alignment and delayed network returns',
                     'Model time is shown; the two solvers execute sequentially.')
    _panel_heading(ax, 'a', 'Fixed exchange intervals', 121,
                   'A terminal interval may be shorter')
    ticks = [255, 485, 715, 945]
    for n, x in enumerate(ticks):
        text(ax, x, 179, rf'$t_{n}$', size=15, ha='center', color=GRAY)
        line(ax, [x, x], [196, 377], GRAY, 0.9, dashed=True)
    text(ax, 35, 222, 'Surface', size=15.5, color=BLUE, weight='bold')
    text(ax, 35, 288, 'SWMM', size=15.5, color=TEAL, weight='bold')
    text(ax, 35, 355, 'Positive\nreturn rate', size=14, color=ORANGE)
    for x, widths in zip(ticks, ([57.5]*4, [230/3]*3, [46]*5)):
        _time_strip(ax, x, 205, widths, BLUE, PALE_BLUE)
    for x, widths in zip(ticks, ([65, 30, 50, 85], [30, 60, 40, 50, 50], [85, 65, 80])):
        _time_strip(ax, x, 271, widths, TEAL, PALE_TEAL)
    rect(ax, 255, 337, 230, 34, '#F2F3F4', edge=LINE, lw=0.9)
    text(ax, 370, 354, 'Initial reservoir', size=13.2, color=GRAY, ha='center')
    for n, x in enumerate((485, 715), start=1):
        rect(ax, x, 337, 230, 34, PALE_ORANGE, edge=ORANGE, lw=1.2)
        text(ax, x+115, 354, rf'$P_s^{n}/T_{n+1}$', size=16, color=ORANGE, ha='center')
        ax.add_patch(FancyArrowPatch((x, 306), (x+34, 334),
                     connectionstyle='arc3,rad=-0.25', arrowstyle='-|>',
                     mutation_scale=14, color=ORANGE, linewidth=1.5, zorder=5))
    text(ax, 255, 408, r'$\sum_j\Delta t_j=\sum_k\delta t_k=T_n$', size=15, color=TEAL)
    text(ax, 965, 408, r'$P_s^n$: queue just after exchange at $t_n$',
         size=13.2, color=GRAY, ha='right')
    text(ax, 255, 451, r'$V_i^n=\sum_j Q_{i,j}\,\Delta t_j$', size=16)
    text(ax, 965, 451, r'$\overline{Q}_i^n=V_i^n/T_n$', size=17, ha='right')
    line(ax, [35, 965], [479, 479], LINE, 0.9)

    _panel_heading(ax, 'b', 'Exchange after every surface step', 513,
                   r'$T_j=\Delta t_j$')
    step_ticks = [255, 445, 685, 945]
    for j, x in enumerate(step_ticks):
        text(ax, x, 552, rf'$t^{{{j}}}$', size=15, ha='center', color=GRAY)
        line(ax, [x, x], [569, 684], GRAY, 0.9, dashed=True)
    text(ax, 35, 593, 'Surface', size=15.5, color=BLUE, weight='bold')
    text(ax, 35, 659, 'SWMM', size=15.5, color=TEAL, weight='bold')
    _time_strip(ax, 255, 576, [190, 240, 260], BLUE, PALE_BLUE)
    for j, (x, width) in enumerate(zip(step_ticks, (190, 240, 260)), start=1):
        text(ax, x+width/2, 593, rf'$\Delta t_{j}$', size=16, color=BLUE, ha='center')
    for x, widths in zip(step_ticks, ([70, 45, 75], [60, 100, 80], [85, 65, 110])):
        _time_strip(ax, x, 642, widths, TEAL, PALE_TEAL)
    for x in step_ticks[1:]:
        arrow(ax, x, 640, x, 612, color=ORANGE, lw=1.6)
    text(ax, 255, 715, 'Orange arrows: return volumes apply in the current surface step.',
         size=13.5, color=GRAY)
    line(ax, [35, 965], [741, 741], LINE, 0.9)

    _panel_heading(ax, 'c', 'Current CFL limit and exact target landing', 777)
    text(ax, 48, 841, r'$N=\left\lceil\frac{R}{\Delta t_{\mathrm{CFL}}}\right\rceil$', size=18)
    text(ax, 48, 894, r'$\Delta t=\frac{R}{N}\leq\Delta t_{\mathrm{CFL}}$', size=19, color=BLUE)
    text(ax, 48, 943, 'Recompute the CFL limit\nas the solution changes.',
         size=13.5, color=GRAY, va='top')
    text(ax, 687, 815, '$R$: time remaining', size=14, color=GRAY, ha='center')
    arrow(ax, 430, 843, 945, 843, color=GRAY, lw=1.2, both=True)
    text(ax, 430, 873, '$t$', size=15, color=GRAY, ha='center')
    text(ax, 945, 873, r'target $t_*$', size=14, color=GRAY, ha='right')
    _time_strip(ax, 430, 895, [103]*5, BLUE, PALE_BLUE)
    arrow(ax, 430, 954, 533, 954, color=BLUE, lw=1.3, both=True)
    text(ax, 481.5, 983, r'$\Delta t$', size=15, color=BLUE, ha='center')
    arrow(ax, 692, 954, 814, 954, color=GRAY, lw=1.3, both=True)
    text(ax, 753, 983, r'$\Delta t_{\mathrm{CFL}}$', size=15, color=GRAY, ha='center')
    text(ax, 35, 1031, 'Targets: exchange boundary, rain frame, requested restart or simulation end.',
         size=13.5, color=GRAY)
    text(ax, 35, 1064, 'Rain-frame alignment tolerance: 0.05 s. Rainfall follows its time coordinates.',
         size=13.2, color=GRAY)
    return fig


FIGURES = {f.__name__: f for f in (surface, layers, exchange_cycle,
                                  manhole_exchange, outfall_cv, time_sync)}


def validate_layout(fig):
    """Reject cropped or overlapping text before exporting the figure."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = fig.axes[0].texts
    for i, label in enumerate(texts):
        box = label.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0, label.get_text()
        assert box.x1 <= fig.bbox.x1 and box.y1 <= fig.bbox.y1, label.get_text()
        anchor = label.get_transform().transform(label.get_position())
        for patch in fig.axes[0].patches:
            if not getattr(patch, '_constrain_text', False):
                continue
            boundary = patch.get_window_extent(renderer)
            if boundary.contains(*anchor):
                assert (box.x0 >= boundary.x0+3 and box.x1 <= boundary.x1-3
                        and box.y0 >= boundary.y0+3 and box.y1 <= boundary.y1-3), label.get_text()
        for other in texts[i+1:]:
            assert not box.overlaps(other.get_window_extent(renderer)), (
                label.get_text(), other.get_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1]/'docs/figures')
    parser.add_argument('--only', choices=list(FIGURES), nargs='+')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in args.only or FIGURES:
        fig = FIGURES[name]()
        validate_layout(fig)
        svg = StringIO()
        fig.savefig(svg, format='svg', metadata={
            'Date': None, 'Title': fig.axes[0].texts[0].get_text(),
            'Description': 'LISFLOOD-FP and SWMM coupling schematic. Source: tools/draw_figures.py.'})
        (args.output/f'{name}.svg').write_text(
            '\n'.join(line.rstrip() for line in svg.getvalue().splitlines()) + '\n',
            encoding='utf-8', newline='\n')
        fig.savefig(args.output/f'{name}.png', dpi=300)
        plt.close(fig)
        print(name, flush=True)


if __name__ == '__main__':
    main()
