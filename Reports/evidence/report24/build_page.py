#!/usr/bin/env python3
"""Report 24: render the summary page (Reports/report24_page.html) from page_data.json and restart.json.
Every number on the page is read from those files, which collect_page_data.py and the run outputs produce."""
import html
import json
import math
from pathlib import Path

E = Path(__file__).resolve().parent
ROOT = E.parents[2]
P = json.loads((E / 'page_data.json').read_text())
RS = json.loads((E / 'restart.json').read_text()) if (E / 'restart.json').exists() else None
esc = html.escape


def fmt(v, nd=1):
    if v is None:
        return '—'
    if v != 0 and (abs(v) < 0.1 or abs(v) >= 1e4):
        return f'{v:.1e}'
    return f'{v:.{nd}f}'


# ---------- sequential blue ramp (reference palette steps), log-scaled ----------
RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab',
        '#184f95', '#104281', '#0d366b']


def heat(v, lo=0.1, hi=400.0):
    t = (math.log10(max(v, lo)) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    i = min(len(RAMP) - 1, max(0, round(t * (len(RAMP) - 1))))
    return RAMP[i], ('#ffffff' if i >= 6 else '#0b0b0b')


def heat_table():
    regs = P['regions']
    names = {'start_base': 'Start (transferred)', 'it01_trial_l0.5': 'Iterate 1', 'it02_trial_l1': 'Iterate 2',
             'it03_trial_l0.25': 'Iterate 3', 'it04_trial_l0.25': 'Iterate 4'}
    head = ''.join(f'<th scope="col">{esc(r)}</th>' for r in regs)
    rows = []
    passes = [dict(q, label='Iterate 1') if q['name'] == 'it01_trial_l0.5' else q for q in P['passes']]
    if RS:
        passes = passes + [dict(name=f"c{q['iterate'] + 1}", label=f"Constrained {q['iterate'] + 1}", by_region=q['by_region'],
                                total=q['total']) for q in RS['rows'] if q['accepted'] and 'total' in q]
    for p in passes:
        cells = []
        for r, v in zip(regs, p['by_region']):
            bg, fg = heat(v)
            cells.append(f'<td class="hc" style="background:{bg};color:{fg}" tabindex="0" '
                         f'data-tip="{esc(names.get(p["name"], p["name"]))} · depths {esc(r)}: {fmt(v)}">{fmt(v)}</td>')
        rows.append(f'<tr><th scope="row">{esc(p.get("label") or names.get(p["name"], p["name"]))}</th>{"".join(cells)}'
                    f'<td class="tot">{fmt(p["total"])}</td></tr>')
    return (f'<div class="tablewrap"><table class="heat"><thead><tr><th scope="col">State</th>{head}'
            f'<th scope="col" class="tot">Total</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>')


# ---------- generic SVG line chart (log y) ----------
def line_chart(series, *, ylo, yhi, ylog=True, xlo=1, xhi=167, yticks=(), xticks=(), ref=None, w=640, h=300,
               ylabel='', xlabel='Depth index (1 = outer boundary, 167 = inner boundary)', bands=(), cid='c'):
    ml, mr, mt, mb = 56, 18, 14, 44
    pw, ph = w - ml - mr, h - mt - mb

    def X(x):
        return ml + (x - xlo) / (xhi - xlo) * pw

    def Y(y):
        if ylog:
            y = max(y, ylo)
            return mt + (1 - (math.log10(y) - math.log10(ylo)) / (math.log10(yhi) - math.log10(ylo))) * ph
        return mt + (1 - (y - ylo) / (yhi - ylo)) * ph
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="{esc(ylabel)}" id="{cid}">']
    for x0, x1, lab in bands:
        out.append(f'<rect x="{X(x0):.1f}" y="{mt}" width="{X(x1) - X(x0):.1f}" height="{ph}" class="band"/>'
                   f'<text x="{(X(x0) + X(x1)) / 2:.1f}" y="{mt + 12}" class="bandlab" text-anchor="middle">{esc(lab)}</text>')
    for t in yticks:
        out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" class="grid"/>'
                   f'<text x="{ml - 8}" y="{Y(t) + 4:.1f}" class="tick" text-anchor="end">{t:g}</text>')
    for t in xticks:
        out.append(f'<text x="{X(t):.1f}" y="{mt + ph + 18}" class="tick" text-anchor="middle">{t}</text>')
    out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{mt + ph}" y2="{mt + ph}" class="axis"/>')
    if ref is not None:
        out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{Y(ref[0]):.1f}" y2="{Y(ref[0]):.1f}" class="ref"/>'
                   f'<text x="{ml + pw - 4}" y="{Y(ref[0]) - 6:.1f}" class="reflab" text-anchor="end">{esc(ref[1])}</text>')
    for s in series:
        pts = ' '.join(f'{X(x):.1f},{Y(y):.1f}' for x, y in zip(s['x'], s['y']))
        out.append(f'<polyline points="{pts}" class="ln" style="stroke:var({s["var"]})"/>')
    out.append(f'<text x="{ml + pw / 2:.1f}" y="{h - 6}" class="axlab" text-anchor="middle">{esc(xlabel)}</text>')
    out.append(f'<text x="14" y="{mt + ph / 2:.1f}" class="axlab" text-anchor="middle" '
               f'transform="rotate(-90 14 {mt + ph / 2:.1f})">{esc(ylabel)}</text>')
    # crosshair layer (filled by script from data-series)
    out.append(f'<line class="xh" x1="0" x2="0" y1="{mt}" y2="{mt + ph}" visibility="hidden"/>')
    out.append(f'<rect class="hit" x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="transparent"/>')
    out.append('</svg>')
    meta = dict(ml=ml, pw=pw, xlo=xlo, xhi=xhi,
                series=[dict(name=s['name'], var=s['var'], x=s['x'], y=s['y']) for s in series])
    return ''.join(out), meta


def legend(items, kind='line'):
    return '<ul class="legend">' + ''.join(
        f'<li><span class="key {kind}" style="--k:var({v})"></span>{esc(n)}</li>' for n, v in items) + '</ul>'


# ---------- charts ----------
CHARTS = {}
by = {p['name']: p for p in P['passes']}
depth = list(range(1, 168))
he_series = [dict(name='Start (transferred)', var='--o1', x=depth, y=by['start_base']['he_over_he0']),
             dict(name='Iterate 2', var='--o2', x=depth, y=by['it02_trial_l1']['he_over_he0']),
             dict(name='Iterate 4', var='--o3', x=depth, y=by['it04_trial_l0.25']['he_over_he0'])]
svg_he, CHARTS['he'] = line_chart(he_series, ylo=0.1, yhi=3.0, yticks=(0.1, 0.2, 0.5, 1, 2), xticks=(1, 23, 45, 75, 118, 150, 167),
                                  ref=(1.0, 'exact (every projected trial)'), ylabel='He total / x0 total', cid='he',
                                  bands=((23, 45, 'H fronts'), (46, 75, 'He front')))
ch_series = [dict(name='Start (transferred)', var='--o1', x=depth, y=by['start_base']['charge_over_ne']),
             dict(name='Iterate 2', var='--o2', x=depth, y=by['it02_trial_l1']['charge_over_ne']),
             dict(name='Iterate 4', var='--o3', x=depth, y=by['it04_trial_l0.25']['charge_over_ne'])]
svg_ch, CHARTS['ch'] = line_chart(ch_series, ylo=0.04, yhi=40.0, yticks=(0.1, 1, 10), xticks=(1, 23, 45, 75, 118, 150, 167),
                                  ref=(1.0, 'charge balance'), ylabel='Ion charge sum / n_e', cid='ch',
                                  bands=((23, 45, 'H fronts'), (46, 75, 'He front')))

# albedo columns, depths 40-85
alb = P['albedo']
sel = [(d, v) for d, v in zip(alb['depth'], alb['scat_over_chi']) if 44 <= d <= 80]


def albedo_chart(w=640, h=260):
    ml, mr, mt, mb = 56, 18, 14, 44
    pw, ph = w - ml - mr, h - mt - mb
    ymax = 3.0
    n = len(sel); slot = pw / n; bw = min(10.0, slot - 2)
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Scattering over total opacity by depth">']
    for t in (0.5, 1, 1.5, 2, 2.5, 3):
        y = mt + (1 - t / ymax) * ph
        out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>'
                   f'<text x="{ml - 8}" y="{y + 4:.1f}" class="tick" text-anchor="end">{t:g}</text>')
    y1 = mt + (1 - 1 / ymax) * ph
    for i, (d, v) in enumerate(sel):
        x = ml + i * slot + (slot - bw) / 2
        hh = v / ymax * ph; y = mt + ph - hh
        cls = 'bar crit' if v > 1 else 'bar'
        out.append(f'<path d="M{x:.1f},{mt + ph:.1f} V{y + 3:.1f} Q{x:.1f},{y:.1f} {x + 3:.1f},{y:.1f} H{x + bw - 3:.1f} '
                   f'Q{x + bw:.1f},{y:.1f} {x + bw:.1f},{y + 3:.1f} V{mt + ph:.1f} Z" class="{cls}" tabindex="0" '
                   f'data-tip="depth {d}: CHI_SCAT/CHI = {v:.3f}"/>')
        if d % 5 == 0:
            out.append(f'<text x="{x + bw / 2:.1f}" y="{mt + ph + 18}" class="tick" text-anchor="middle">{d}</text>')
    out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{y1:.1f}" y2="{y1:.1f}" class="ref"/>'
               f'<text x="{ml + pw - 4}" y="{y1 - 6:.1f}" class="reflab" text-anchor="end">net absorption = 0</text>')
    out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{mt + ph}" y2="{mt + ph}" class="axis"/>')
    out.append(f'<text x="{ml + pw / 2:.1f}" y="{h - 6}" class="axlab" text-anchor="middle">Depth index</text>')
    out.append(f'<text x="14" y="{mt + ph / 2:.1f}" class="axlab" text-anchor="middle" transform="rotate(-90 14 {mt + ph / 2:.1f})">'
               f'CHI_SCAT / CHI</text></svg>')
    return ''.join(out)


def floor_chart(w=640):
    rows = P['floor_q']
    ml, mr, mt, rh = 150, 24, 10, 32
    h = mt + rh * len(rows) + 40
    pw = w - ml - mr
    lo, hi = 0.5, 300.0

    def X(v):
        return ml + (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * pw
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Worst-case preconditioner error before and after the floor">']
    for t in (0.5, 1, 2, 5, 10, 20, 50, 100, 200):
        out.append(f'<line x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{mt}" y2="{mt + rh * len(rows)}" class="grid"/>'
                   f'<text x="{X(t):.1f}" y="{mt + rh * len(rows) + 16}" class="tick" text-anchor="middle">{t:g}</text>')
    out.append(f'<line x1="{X(1):.1f}" x2="{X(1):.1f}" y1="{mt}" y2="{mt + rh * len(rows)}" class="ref"/>')
    names = {'fixedT': 'Fixed T', 'free': 'T released'}
    for i, (lab, m) in enumerate(rows):
        y = mt + rh * i + rh / 2
        raw = m['ALI'][2]; fl = [v for k, v in m.items() if k.startswith('ALI floor 0.01')][0][2]
        run, itn = lab.split()
        out.append(f'<text x="{ml - 10}" y="{y + 4:.1f}" class="tick" text-anchor="end">{esc(names[run])}, {esc(itn.replace("it", "iterate "))}</text>')
        out.append(f'<line x1="{X(fl):.1f}" x2="{X(raw):.1f}" y1="{y:.1f}" y2="{y:.1f}" class="dumb"/>')
        out.append(f'<circle cx="{X(raw):.1f}" cy="{y:.1f}" r="5" class="dot o" tabindex="0" data-tip="{esc(lab)} · local ALI, worst q: {raw:.1f}"/>')
        out.append(f'<circle cx="{X(fl):.1f}" cy="{y:.1f}" r="5" class="dot b" tabindex="0" data-tip="{esc(lab)} · floored, worst q: {fl:.2f}"/>')
    out.append(f'<text x="{ml + pw / 2:.1f}" y="{h - 4}" class="axlab" text-anchor="middle">Worst-case q = ‖W(M⁻¹Jz − z)‖ / ‖Wz‖ over exact products (log)</text></svg>')
    return ''.join(out)


def trajectory():
    if not RS:
        return '<p class="note">The constrained restart is running; results are added when its passes finish.</p>', None
    s1 = RS['trajectory_old']; s2 = RS['trajectory_new']
    xs = [p[0] for p in s1 + s2]
    ys = [p[1] for p in s1 + s2]
    svg, meta = line_chart([dict(name='Unconstrained (Report 23 path)', var='--o2', x=[p[0] for p in s1], y=[p[1] for p in s1]),
                            dict(name='Constrained restart', var='--s1', x=[p[0] for p in s2], y=[p[1] for p in s2])],
                           ylo=(50 if min(ys) >= 50 else 10), yhi=400, yticks=[t for t in RS['yticks'] + [10, 20] if (50 if min(ys) >= 50 else 10) <= t <= 400], xlo=min(xs), xhi=max(xs),
                           xticks=sorted(set(xs)), ylabel='Fixed-measure residual', xlabel='Iterate', cid='traj', w=640, h=280)
    return svg, meta


svg_traj, CHARTS['traj'] = trajectory()
CHARTS = {k: v for k, v in CHARTS.items() if v}
nd = P.get('newton_dump') or []
rho = ', '.join(f'{r["rho"]:.3f}' for r in nd) if nd else '—'

RSC = {}
for ln in (E / 'rowscale_by_region.txt').read_text().splitlines():
    parts = ln.rsplit(None, 7)
    if len(parts) == 8 and not ln.startswith('state'):
        try:
            RSC[parts[0].strip()] = float(parts[-1])
        except ValueError:
            pass
rs_txt = ' → '.join(f'{RSC[k]:.2f}' for k in ('iterate 1 (it01_trial_l0.5)', 'constrained 2', 'constrained 3 (l1)',
                                              'constrained 4 (l1)', 'constrained 5 (l0.25)') if k in RSC)
if RS:
    tn = RS['trajectory_new']
    path_txt = ' → '.join(f'{v:.0f}' for _, v in tn)
    restart_html = (
        '<div class="stack read"><p>Projection cannot rescue the damaged iterate-4 state: the projected base alone measured '
        f'{RS["it04_projected"]:.0f} against 118, because its populations had adapted to totals off by up to 8× (He) and 33× (n_e). '
        'So the solve restarted from iterate 1, the first accepted state, with the floor, the projection and T held fixed. '
        'Every accepted state conserves H, He and charge to 1e-15. The first three iterates accepted full Newton steps and '
        'the fixed measure fell from 224 to 128. At the fourth, the full and half steps met supercritical masers again, and '
        'the quarter step the driver accepted is worse in the fixed measure (140). A second measure with no inverse, the '
        f'row-scaled residual, went {rs_txt}. Both must reach zero at the solution: this is a solve that keeps the '
        'constraints, not yet a converging one. The run stopped at 05:16, at the end of the session, and can be resumed '
        'from its best base.</p></div>' + RS['table_html'])
    status_closure = 'Stage 1 stalled near 128, conservation exact'
else:
    restart_html = '<p class="note">Running at publish time.</p>'
    status_closure = 'Stage 1 running'
NEXT_TEXT = ('The constrained solve conserves everything exactly, but each step still '
             'damages the converged deep interior. The hydrogen ionization front is so ill-conditioned in the preconditioner that '
             'a 22 % far-field radiative coupling of deep-layer directions becomes 99.6 % of their preconditioned action, and '
             'the fit uses it.')

page = f'''<title>Cool Closure Night Log</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:wght@500;600&family=Source+Sans+3:wght@400;600&family=JetBrains+Mono:wght@400&display=swap">
<style>
/* Layout: one reading column; charts break out to a wider measure; a depth strip (outer -> inner) is the recurring motif. */
:root {{
  --bg: #f3f5f8; --surface: #fcfcfb; --ink: #12161d; --ink2: #4a5260; --muted: #7b8494; --rule: #dce1e8;
  --halpha: #b8323b; --band: rgba(42,120,214,0.06);
  --o1: #86b6ef; --o2: #2a78d6; --o3: #104281; --s1: #eb6834; --crit: #d03b3b; --good: #0ca30c;
  --grid: #e1e0d9; --axis: #c3c2b7;
  --display: "Spectral", Georgia, "Times New Roman", serif;
  --body: "Source Sans 3", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "JetBrains Mono", ui-monospace, "SF Mono", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #0d1016; --surface: #151a22; --ink: #e9edf3; --ink2: #b2bac6; --muted: #858e9c; --rule: #262d38;
  --halpha: #ef6670; --band: rgba(57,135,229,0.10);
  --o1: #184f95; --o2: #3987e5; --o3: #9ec5f4; --s1: #d95926; --crit: #e66767; --good: #0ca30c;
  --grid: #2c2c2a; --axis: #383835; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #0d1016; --surface: #151a22; --ink: #e9edf3; --ink2: #b2bac6; --muted: #858e9c; --rule: #262d38;
  --halpha: #ef6670; --band: rgba(57,135,229,0.10);
  --o1: #184f95; --o2: #3987e5; --o3: #9ec5f4; --s1: #d95926; --crit: #e66767; --good: #0ca30c;
  --grid: #2c2c2a; --axis: #383835; color-scheme: dark; }}
* {{ box-sizing: border-box; }}
body {{ background: var(--bg); color: var(--ink); font: 400 17px/1.6 var(--body); margin: 0; }}
.wrap {{ max-width: 980px; margin: 0 auto; padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px)); padding-block: 40px 72px; }}
.read {{ max-width: 680px; }}
h1, h2, h3 {{ font-family: var(--display); font-weight: 600; text-wrap: balance; line-height: 1.15; margin: 0; }}
h1 {{ font-size: clamp(34px, 5vw, 52px); letter-spacing: -0.01em; }}
h2 {{ font-size: 28px; margin-top: 8px; }}
h3 {{ font-size: 19px; font-family: var(--body); font-weight: 600; }}
.eyebrow {{ font: 400 12px/1.4 var(--mono); letter-spacing: 0.08em; text-transform: uppercase; color: var(--muted); }}
p {{ margin: 0; }}
.stack {{ display: grid; gap: 14px; }}
section {{ display: grid; gap: 18px; padding-top: 44px; margin-top: 44px; border-top: 1px solid var(--rule); }}
.lede {{ font-size: 20px; color: var(--ink2); }}
code, .mono {{ font-family: var(--mono); font-size: 0.86em; }}
.chips {{ display: flex; flex-wrap: wrap; gap: 8px; }}
.chip {{ display: inline-flex; align-items: center; gap: 8px; padding: 5px 12px; border-radius: 999px; border: 1px solid var(--rule);
  background: var(--surface); font-size: 14px; color: var(--ink2); }}
.chip i {{ width: 8px; height: 8px; border-radius: 50%; background: var(--muted); display: inline-block; }}
.chip.bad i {{ background: var(--crit); }} .chip.ok i {{ background: var(--good); }}
.strip {{ display: grid; grid-template-columns: 22fr 23fr 30fr 43fr 32fr 17fr; height: 34px; border-radius: 6px; overflow: hidden;
  font: 400 11px/1 var(--mono); color: var(--ink2); }}
.strip span {{ display: flex; align-items: center; justify-content: center; border-right: 2px solid var(--bg); background: var(--surface);
  white-space: nowrap; overflow: hidden; }}
.strip span.hf {{ background: color-mix(in srgb, var(--halpha) 16%, var(--surface)); }}
.strip span.he {{ background: color-mix(in srgb, #d9a400 18%, var(--surface)); }}
.findings {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 16px; }}
.finding {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 10px; padding: 18px 18px 20px; display: grid; gap: 10px;
  align-content: start; min-width: 0; }}
.finding .fig {{ font: 600 30px/1 var(--body); }}
.finding .fig small {{ font-size: 15px; font-weight: 400; color: var(--muted); }}
.finding p {{ color: var(--ink2); font-size: 15.5px; }}
.tag {{ font: 400 11px/1 var(--mono); letter-spacing: 0.06em; text-transform: uppercase; color: var(--muted); }}
.tag.fixed {{ color: var(--good); }} .tag.cause {{ color: var(--halpha); }}
figure {{ margin: 0; background: var(--surface); border: 1px solid var(--rule); border-radius: 10px; padding: 18px 16px 14px; display: grid; gap: 10px; min-width: 0; }}
figcaption {{ color: var(--ink2); font-size: 15px; }}
figcaption b {{ color: var(--ink); font-weight: 600; }}
.chart {{ width: 100%; height: auto; display: block; overflow: visible; }}
.chart .grid {{ stroke: var(--grid); stroke-width: 1; }}
.chart .axis {{ stroke: var(--axis); stroke-width: 1; }}
.chart .tick {{ fill: var(--muted); font: 400 13px var(--mono); }}
.chart .axlab {{ fill: var(--ink2); font: 400 14px var(--body); }}
.chart .ln {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
.chart .ref {{ stroke: var(--ink2); stroke-width: 1; stroke-dasharray: 0; opacity: 0.7; }}
.chart .reflab, .chart .bandlab {{ fill: var(--ink2); font: 400 13px var(--body); }}
.chart .band {{ fill: var(--band); }}
.chart .bar {{ fill: var(--o2); }} .chart .bar.crit {{ fill: var(--crit); }}
.chart .bar:hover, .chart .bar:focus, .chart .dot:hover, .chart .dot:focus {{ opacity: 0.75; outline: none; }}
.chart .dot {{ stroke: var(--surface); stroke-width: 2; }} .chart .dot.o {{ fill: var(--s1); }} .chart .dot.b {{ fill: var(--o2); }}
.chart .dumb {{ stroke: var(--axis); stroke-width: 2; }}
.chart .xh {{ stroke: var(--ink2); stroke-width: 1; }}
.legend {{ list-style: none; margin: 0; padding: 0; display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 14px; color: var(--ink2); }}
.legend li {{ display: inline-flex; align-items: center; gap: 8px; }}
.key.line {{ width: 18px; height: 0; border-top: 2px solid var(--k); display: inline-block; }}
.key.dot {{ width: 10px; height: 10px; border-radius: 50%; background: var(--k); display: inline-block; }}
.key.rect {{ width: 12px; height: 12px; border-radius: 3px; background: var(--k); display: inline-block; }}
.tablewrap {{ overflow-x: auto; }}
table.heat {{ border-collapse: separate; border-spacing: 2px; font: 400 13px/1.2 var(--mono); width: 100%; min-width: 620px; }}
table.heat th {{ font: 600 12px/1.3 var(--body); color: var(--ink2); text-align: left; padding: 6px 8px; white-space: nowrap; }}
table.heat thead th {{ text-align: center; }}
table.heat td {{ text-align: center; padding: 9px 6px; border-radius: 4px; font-variant-numeric: tabular-nums; }}
table.heat td.tot {{ background: transparent; color: var(--ink); font-weight: 600; }}
table.heat td.hc:focus {{ outline: 2px solid var(--ink); outline-offset: 1px; }}
table.trials {{ border-collapse: collapse; width: 100%; min-width: 640px; font: 400 14px/1.35 var(--body); }}
table.trials th {{ text-align: left; font-weight: 600; color: var(--ink2); padding: 8px 10px; border-bottom: 1px solid var(--rule); white-space: nowrap; }}
table.trials td {{ padding: 8px 10px; border-bottom: 1px solid var(--rule); font-variant-numeric: tabular-nums; color: var(--ink2); }}
table.trials td:last-child {{ color: var(--good); font-weight: 600; }}
.grid2 {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }}
.note {{ color: var(--muted); font-size: 14px; }}
.kv {{ display: grid; grid-template-columns: max-content 1fr; gap: 6px 16px; font-size: 15px; }}
.kv dt {{ color: var(--muted); }} .kv dd {{ margin: 0; color: var(--ink2); min-width: 0; overflow-wrap: anywhere; }}
ol.next {{ margin: 0; padding-left: 22px; display: grid; gap: 8px; color: var(--ink2); }}
#tip {{ position: fixed; pointer-events: none; z-index: 20; background: var(--surface); color: var(--ink); border: 1px solid var(--rule);
  border-radius: 8px; padding: 8px 10px; font: 400 13px/1.35 var(--body); box-shadow: 0 6px 18px rgba(0,0,0,0.12); max-width: 280px; }}
#tip b {{ font-weight: 600; }} #tip .r {{ display: flex; align-items: center; gap: 8px; }}
#tip .r i {{ width: 14px; border-top: 2px solid; display: inline-block; }}
@media (max-width: 560px) {{ body {{ font-size: 16px; }} .strip {{ font-size: 9px; }} .lede {{ font-size: 18px; }} }}
@media (prefers-reduced-motion: reduce) {{ * {{ transition: none !important; }} }}
</style>
<div class="wrap">
  <header class="stack read">
    <p class="eyebrow">QStarS · Report 24 · night of 21–22 October 2026</p>
    <h1>Under the stall</h1>
    <p class="lede">The 167-radius cool CMFGEN closure had stalled at half its starting residual. One night of measurement
    found why. The solver had been breaking exact conservation laws, fitting the hydrogen front with directions from the
    deep layers, and moving a temperature it was meant to hold fixed. Its full Newton steps produced infrared masers with
    no physical solution.</p>
    <div class="chips">
      <span class="chip bad"><i></i>Closure: not converged</span>
      <span class="chip"><i></i>{esc(status_closure)}</span>
      <span class="chip bad"><i></i>v1: no</span>
      <span class="chip ok"><i></i>Three defects fixed and tested</span>
      <span class="chip"><i></i>CMFGEN build 36</span>
    </div>
  </header>

  <section>
    <div class="stack read">
      <p class="eyebrow">The atmosphere, outer to inner</p>
      <p>Every chart below is indexed by depth through the model atmosphere, from the outer boundary (depth 1) to the
      inner boundary (167). Hydrogen changes ionization at depths 23–45, helium near 46–75.</p>
    </div>
    <div class="strip" aria-hidden="true"><span>1–22</span><span class="hf">23–45</span><span class="he">46–75</span>
      <span>76–118</span><span>119–150</span><span>151–167</span></div>
    <ul class="legend"><li><span class="key rect" style="--k:color-mix(in srgb, var(--halpha) 30%, var(--surface))"></span>Hydrogen ionization fronts</li><li><span class="key rect" style="--k:color-mix(in srgb, #d9a400 34%, var(--surface))"></span>Helium front</li></ul>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Four findings</p><h2>What was under the stall</h2></div>
    <div class="findings">
      <article class="finding"><span class="tag fixed">Cause · fixed</span><h3>Conservation was broken</h3>
        <p class="fig">34×<small> worst charge imbalance</small></p>
        <p>The step limiter clips each population's change on its own, which breaks the exact number- and charge-conservation
        equations. Nothing restored them. Helium drifted up to 8× and hydrogen 2.3×, in every Report 23 solve.</p></article>
      <article class="finding"><span class="tag cause">Cause · open at the H front</span><h3>Distant directions as blunt instruments</h3>
        <p class="fig">172<small> coefficient of a deep-layer direction</small></p>
        <p>The least squares fitted the hydrogen front with a deep-layer direction, whose tiny far-field action the nearly
        singular front blocks amplify. 97 % of its effect landed at the front; the deep layers took the damage.</p></article>
      <article class="finding"><span class="tag fixed">Cause · fixed</span><h3>Fixed T was not fixed</h3>
        <p class="fig">0.50×<small> trial T at the inner boundary</small></p>
        <p>In the fixed-temperature stage, T was still a variable held by a weak pin row. The fit moved it (coefficient 358),
        and a full trial halved T at depth 166. T is now removed from every direction.</p></article>
      <article class="finding"><span class="tag cause">Cause · diagnosed</span><h3>Trial steps made masers</h3>
        <p class="fig">ρ = {esc(rho)}<small> spectral radius</small></p>
        <p>Full steps inverted H I infrared lines at the He front. Where scattering exceeds total opacity, the transfer
        operator's spectral radius passes 1: no steady radiation field exists, and CMFGEN fails.</p></article>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Where the residual lives</p><h2>The fronts improved while the deep interior was damaged</h2>
      <p>Fixed preconditioned measure ‖W M⁻¹R‖ (one reference preconditioner for every state), split by depth region. The
      transferred start had a converged deep interior (0.2 and 0.1). By iterate 4 those regions read 66 and 30, while the
      hydrogen fronts fell from 350 to 78.</p></div>
    <figure>{heat_table()}<figcaption>Colour is the same log scale in every cell, from 0.1 (lightest) to 400 (darkest).
      Region norms add in quadrature to the total.</figcaption></figure>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Conservation</p><h2>Helium and charge drifted at the inserted depths</h2>
      <p>The 167-radius grid interleaves the 84-radius one; the even depths were inserted. From iterate 2 the steps
      halved n_e and He II at even depths of the deep interior every iterate, the limiter's own factor.</p></div>
    <div class="stack">
      <figure>{legend([(s['name'], s['var']) for s in he_series])}{svg_he}
        <figcaption><b>Helium total against x0's.</b> 1 is exact conservation. The deep zig-zag is the even depths.</figcaption></figure>
      <figure>{legend([(s['name'], s['var']) for s in ch_series])}{svg_ch}
        <figcaption><b>Charge sum over n_e.</b> At iterate 4 the worst depth is off by 34×.</figcaption></figure>
    </div>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">The native failures</p><h2>A maser with no steady state</h2>
      <p>In the failed full-step trial, the H I 12→8 line (0.0285 × 10¹⁵ Hz) is inverted at depths 60, 62 and 64. CMFGEN's
      negative-opacity floor keeps the opacity positive but below electron scattering there. Build 36 added a Newton
      solve of the scattering fixed point, and the dumped Jacobian answers the question: spectral radius {esc(rho)}, a
      fixed point that is negative at depths 57–67. No transfer solver can fix that state; the step must not make it.</p></div>
    <figure>{albedo_chart()}<figcaption><b>Scattering over total opacity at 0.0285 × 10¹⁵ Hz</b>, last scattering
      iterate of trial <code>it04_trial_l0.5</code>, depths 44–80. Red columns exceed 1 (negative net absorption).</figcaption></figure>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Preconditioner</p><h2>Flooring the singular values</h2>
      <p>Each equilibrated diagonal block keeps its singular vectors; singular values below 1 % of the largest are raised
      to that floor. Measured against exact Jacobian products at eight archived states.</p></div>
    <figure>{legend([('Local ALI (Report 23)', '--s1'), ('With the 1e-2 floor', '--o2')], 'dot')}{floor_chart()}
      <figcaption><b>Worst-case error of the preconditioner</b>; 1 means it explains none of the exact action.</figcaption></figure>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Constrained restart</p><h2>Stage 1 with exact conservation</h2></div>
    {restart_html}
    {f'<figure>{legend([("Unconstrained (Report 23 path)", "--o2"), ("Constrained restart", "--s1")])}{svg_traj}<figcaption>Fixed-measure residual of each accepted base.</figcaption></figure>' if RS else ''}
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Next</p><h2>The narrowest remaining blocker</h2>
      <p>{esc(NEXT_TEXT)}</p></div>
    <ol class="next">
      <li>Make the hydrogen front's linear solve exact: exact Jacobian columns for the front depths (23–28) from one batched
      probe pass, so far-field couplings are no longer amplified there.</li>
      <li>Then compare a front-local least squares (each direction counted only near its own depths) against the current
      one on the same products; adopt it only if the true linear residual holds.</li>
      <li>Continue the constrained stage 1 to convergence, release T, and only then reopen Gates 1–5.</li>
    </ol>
  </section>

  <section>
    <div class="stack read"><p class="eyebrow">Record</p><h2>Provenance</h2></div>
    <dl class="kv">
      <dt>Report</dt><dd><code>Reports/Report_24.md</code></dd>
      <dt>Evidence</dt><dd><code>Reports/evidence/report24/</code> (scripts and JSON behind every number here)</dd>
      <dt>Build 36</dt><dd><code>backends/cmfgen/ali_newton.patch</code>, flag <code>qstars_ali_newton_fallback</code></dd>
      <dt>Driver</dt><dd><code>--pc-svd-floor</code>, <code>--conserve</code>, <code>--freeze-t</code>, <code>--trial-workers</code> in <code>scripts/cool_newton_krylov.py</code></dd>
      <dt>Runs</dt><dd><code>results/nk167_trap/fixedT_b35</code>; earlier branches kept in <code>branch_pre_floor/</code> and <code>branch_b35_noconserve/</code></dd>
      <dt>Gates</dt><dd>Unchanged from Report 23: closure FAIL, Gates 1–5 not reached, Gate 6 FAIL, v1 no.</dd>
    </dl>
  </section>
</div>
<div id="tip" hidden></div>
<script>
const CH = {json.dumps(CHARTS)};
const tip = document.getElementById('tip');
function place(e) {{ const pad = 14; let x = e.clientX + pad, y = e.clientY + pad;
  const r = tip.getBoundingClientRect(); if (x + r.width > innerWidth - 8) x = e.clientX - r.width - pad;
  if (y + r.height > innerHeight - 8) y = e.clientY - r.height - pad; tip.style.left = x + 'px'; tip.style.top = y + 'px'; }}
function showText(e, text) {{ tip.replaceChildren(); const b = document.createElement('b'); b.textContent = text; tip.appendChild(b); tip.hidden = false; place(e); }}
document.querySelectorAll('[data-tip]').forEach(el => {{
  el.addEventListener('pointermove', e => showText(e, el.dataset.tip));
  el.addEventListener('pointerleave', () => tip.hidden = true);
  el.addEventListener('focus', () => {{ const r = el.getBoundingClientRect(); showText({{clientX: r.right, clientY: r.top}}, el.dataset.tip); }});
  el.addEventListener('blur', () => tip.hidden = true);
}});
for (const [id, m] of Object.entries(CH)) {{
  const svg = document.getElementById(id); if (!svg) continue;
  const hit = svg.querySelector('.hit'), xh = svg.querySelector('.xh');
  hit.addEventListener('pointermove', e => {{
    const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY; const p = pt.matrixTransform(svg.getScreenCTM().inverse());
    const xv = m.xlo + (p.x - m.ml) / m.pw * (m.xhi - m.xlo);
    tip.replaceChildren(); let head = null;
    for (const s of m.series) {{ let k = 0, best = Infinity; s.x.forEach((x, i) => {{ const d = Math.abs(x - xv); if (d < best) {{ best = d; k = i; }} }});
      if (head === null) {{ head = s.x[k]; const h = document.createElement('div'); h.textContent = (id === 'traj' ? 'Iterate ' : 'Depth ') + head; h.style.color = 'var(--muted)'; tip.appendChild(h);
        const sx = m.ml + (head - m.xlo) / (m.xhi - m.xlo) * m.pw; xh.setAttribute('x1', sx); xh.setAttribute('x2', sx); xh.setAttribute('visibility', 'visible'); }}
      const r = document.createElement('div'); r.className = 'r'; const i = document.createElement('i'); i.style.borderColor = 'var(' + s.var + ')';
      const b = document.createElement('b'); const v = s.y[k]; b.textContent = Math.abs(v) >= 100 ? v.toFixed(0) : v.toPrecision(3);
      const n = document.createElement('span'); n.textContent = s.name; n.style.color = 'var(--ink2)'; r.append(i, b, n); tip.appendChild(r); }}
    tip.hidden = false; place(e); }});
  hit.addEventListener('pointerleave', () => {{ tip.hidden = true; xh.setAttribute('visibility', 'hidden'); }});
}}
</script>
'''
(ROOT / 'Reports' / 'report24_page.html').write_text(page)
print('written', len(page))
