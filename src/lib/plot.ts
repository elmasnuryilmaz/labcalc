// Publication-style SVG figures, built as strings so the same markup is shown
// in the page and exported as SVG/PNG. Colours: validated categorical palette
// (fixed order, never cycled); text always in neutral ink.
import { logistic4, type FitResult } from './fit4pl';
import { mean, sd, sem, stars } from './stats';

export const PALETTE = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];
const INK = '#0b0b0b';
const INK2 = '#52514e';
const GRID = '#e4e3df';
const FONT = 'Arial, Helvetica, sans-serif';

export const colorFor = (i: number) => (i < PALETTE.length ? PALETTE[i] : '#8a8984');

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const f = (v: number) => (Math.round(v * 100) / 100).toString();

function niceTicks(lo: number, hi: number, target = 6): number[] {
  if (!(hi > lo)) hi = lo + 1;
  const raw = (hi - lo) / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) ?? 10 * mag;
  const out: number[] = [];
  const start = Math.floor(lo / step + 1e-9) * step;
  const end = Math.ceil(hi / step - 1e-9) * step;
  for (let v = start; v <= end + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}

function fmtTick(v: number): string {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1e4 || a < 1e-3) return v.toExponential(0);
  return Number(v.toPrecision(6)).toString();
}

function svgOpen(w: number, h: number, title: string): string {
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="${FONT}" role="img" aria-label="${esc(title || 'figure')}">\n<rect width="${w}" height="${h}" fill="#ffffff"/>\n`;
}

export interface BarGroup {
  group: string;
  values: number[]; // already on the plotted scale
  p?: number; // p-value vs control, for stars
}
export interface BarSeries {
  gene: string;
  groups: BarGroup[];
}
export interface BarOptions {
  log2: boolean;
  error: 'sd' | 'sem';
  showStars: boolean;
  showPoints: boolean;
  title: string;
  yLabel: string;
}

/** Grouped bar chart (gene clusters, one bar per group) with individual biological replicates. */
export function barPlot(series: BarSeries[], groupOrder: string[], o: BarOptions): string {
  const barW = 34;
  const gap = 6;
  const clusterGap = 34;
  const nG = groupOrder.length;
  const clusterW = nG * barW + (nG - 1) * gap;
  const m = { l: 70, r: 20, t: o.title ? 58 : 36, b: 92 };
  const plotW = Math.max(220, series.length * clusterW + (series.length + 1) * clusterGap);
  const legendW = 150;
  const W = m.l + plotW + m.r + legendW;
  const H = 400;
  const plotH = H - m.t - m.b;

  // y range from bars, error bars and points
  let lo = 0;
  let hi = o.log2 ? 0 : 1;
  for (const s of series)
    for (const g of s.groups) {
      if (!g.values.length) continue;
      const mu = mean(g.values);
      const e = (o.error === 'sd' ? sd(g.values) : sem(g.values)) || 0;
      lo = Math.min(lo, mu - e, ...(o.showPoints ? g.values : []));
      hi = Math.max(hi, mu + e, ...(o.showPoints ? g.values : []));
    }
  const pad = (hi - lo) * (o.showStars ? 0.16 : 0.06);
  hi += pad;
  if (lo < 0) lo -= (hi - lo) * 0.04;
  const ticks = niceTicks(lo, hi);
  const y0 = Math.min(lo, ticks[0]);
  const y1 = Math.max(hi, ticks[ticks.length - 1]);
  const Y = (v: number) => m.t + plotH - ((v - y0) / (y1 - y0)) * plotH;

  let s = svgOpen(W, H, o.title);
  if (o.title) s += `<text x="${m.l}" y="24" font-size="15" font-weight="bold" fill="${INK}">${esc(o.title)}</text>\n`;
  for (const t of ticks) {
    s += `<line x1="${m.l}" x2="${m.l + plotW}" y1="${f(Y(t))}" y2="${f(Y(t))}" stroke="${GRID}" stroke-width="1"/>`;
    s += `<text x="${m.l - 8}" y="${f(Y(t) + 4)}" font-size="12" text-anchor="end" fill="${INK2}">${fmtTick(t)}</text>\n`;
  }
  const baseVal = o.log2 ? 0 : 1;
  s += `<line x1="${m.l}" x2="${m.l + plotW}" y1="${f(Y(baseVal))}" y2="${f(Y(baseVal))}" stroke="${INK2}" stroke-width="1" stroke-dasharray="4 3"/>\n`;
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + plotH}" stroke="${INK}" stroke-width="1"/>`;
  s += `<line x1="${m.l}" x2="${m.l + plotW}" y1="${f(Y(Math.max(y0, Math.min(0, y1))))}" y2="${f(Y(Math.max(y0, Math.min(0, y1))))}" stroke="${INK}" stroke-width="1"/>\n`;
  s += `<text transform="translate(18 ${m.t + plotH / 2}) rotate(-90)" font-size="13" text-anchor="middle" fill="${INK}">${esc(o.yLabel)}</text>\n`;

  series.forEach((ser, si) => {
    const cx0 = m.l + clusterGap + si * (clusterW + clusterGap);
    s += `<text x="${f(cx0 + clusterW / 2)}" y="${H - 18}" font-size="13" font-style="italic" text-anchor="middle" fill="${INK}">${esc(ser.gene)}</text>\n`;
    groupOrder.forEach((grp, gi) => {
      const g = ser.groups.find((x) => x.group === grp);
      const x = cx0 + gi * (barW + gap);
      const lx = x + barW / 2;
      s += `<text transform="translate(${f(lx + 4)} ${m.t + plotH + 12}) rotate(-40)" font-size="11" text-anchor="end" fill="${INK2}">${esc(grp)}</text>`;
      if (!g || !g.values.length) return;
      const mu = mean(g.values);
      const e = (o.error === 'sd' ? sd(g.values) : sem(g.values)) || 0;
      const zero = Y(Math.max(y0, Math.min(0, y1)));
      const top = Y(mu);
      const yTop = Math.min(top, zero);
      const hgt = Math.abs(zero - top);
      const tip = `${ser.gene}, ${grp}: mean ${mu.toFixed(3)}, ${o.error.toUpperCase()} ${Number.isFinite(e) ? e.toFixed(3) : 'n/a'}, n = ${g.values.length}`;
      s += `<g><title>${esc(tip)}</title><rect x="${f(x)}" y="${f(yTop)}" width="${barW}" height="${f(Math.max(hgt, 0.5))}" fill="${colorFor(gi)}" fill-opacity="0.85"/>`;
      if (e > 0) {
        s += `<line x1="${f(lx)}" x2="${f(lx)}" y1="${f(Y(mu - e))}" y2="${f(Y(mu + e))}" stroke="${INK}" stroke-width="1.5"/>`;
        s += `<line x1="${f(lx - 7)}" x2="${f(lx + 7)}" y1="${f(Y(mu + e))}" y2="${f(Y(mu + e))}" stroke="${INK}" stroke-width="1.5"/>`;
        s += `<line x1="${f(lx - 7)}" x2="${f(lx + 7)}" y1="${f(Y(mu - e))}" y2="${f(Y(mu - e))}" stroke="${INK}" stroke-width="1.5"/>`;
      }
      if (o.showPoints)
        g.values.forEach((v, k) => {
          const jitter = g.values.length > 1 ? ((k / (g.values.length - 1)) - 0.5) * (barW * 0.5) : 0;
          s += `<circle cx="${f(lx + jitter)}" cy="${f(Y(v))}" r="4" fill="${INK}" fill-opacity="0.8" stroke="#ffffff" stroke-width="1.2"/>`;
        });
      s += `</g>\n`;
      if (o.showStars && g.p !== undefined && Number.isFinite(g.p)) {
        const peak = Math.max(mu + e, mu, ...(o.showPoints ? g.values : []));
        s += `<text x="${f(lx)}" y="${f(Y(peak) - 8)}" font-size="13" text-anchor="middle" fill="${INK}">${stars(g.p)}</text>\n`;
      }
    });
  });

  // legend
  const lx0 = m.l + plotW + m.r + 8;
  groupOrder.forEach((grp, gi) => {
    const ly = m.t + 6 + gi * 22;
    s += `<rect x="${lx0}" y="${ly}" width="14" height="14" rx="2" fill="${colorFor(gi)}" fill-opacity="0.85"/><text x="${lx0 + 20}" y="${ly + 12}" font-size="12" fill="${INK}">${esc(grp)}</text>\n`;
  });
  const noteY = m.t + 6 + groupOrder.length * 22 + 14;
  s += `<text x="${lx0}" y="${noteY}" font-size="11" fill="${INK2}">Bars: mean ± ${o.error.toUpperCase()}</text>`;
  if (o.showPoints) s += `<text x="${lx0}" y="${noteY + 15}" font-size="11" fill="${INK2}">Dots: biological replicates</text>`;
  if (o.showStars) s += `<text x="${lx0}" y="${noteY + 30}" font-size="11" fill="${INK2}">Stars: vs control</text>`;
  return s + '</svg>\n';
}

export interface CurveSeries {
  name: string;
  points: { conc: number; values: number[]; mean: number; sd: number }[];
  fit: FitResult | null;
}
export interface CurveOptions {
  unit: string;
  individual: boolean;
  title: string;
}

/** Dose-response plot on a log10 concentration axis with fitted 4PL curves. */
export function doseResponsePlot(series: CurveSeries[], o: CurveOptions): string {
  const m = { l: 70, r: 20, t: o.title ? 58 : 30, b: 60 };
  const plotW = 460;
  const plotH = 300;
  const legendW = 230;
  const W = m.l + plotW + m.r + legendW;
  const H = m.t + plotH + m.b;
  const allC = series.flatMap((s) => s.points.map((p) => p.conc)).filter((c) => c > 0);
  if (allC.length === 0) return svgOpen(W, H, o.title) + `<text x="20" y="40" font-size="13">No data to plot.</text></svg>`;
  const lx0 = Math.log10(Math.min(...allC)) - 0.35;
  const lx1 = Math.log10(Math.max(...allC)) + 0.35;
  const allY = series.flatMap((s) => s.points.flatMap((p) => (o.individual ? p.values : [p.mean - (p.sd || 0), p.mean + (p.sd || 0)])));
  const ticks = niceTicks(Math.min(0, ...allY), Math.max(110, ...allY), 6);
  const y0 = ticks[0];
  const y1 = ticks[ticks.length - 1];
  const X = (lc: number) => m.l + ((lc - lx0) / (lx1 - lx0)) * plotW;
  const Y = (v: number) => m.t + plotH - ((v - y0) / (y1 - y0)) * plotH;

  let s = svgOpen(W, H, o.title);
  if (o.title) s += `<text x="${m.l}" y="24" font-size="15" font-weight="bold" fill="${INK}">${esc(o.title)}</text>\n`;
  s += `<defs><clipPath id="plotclip"><rect x="${m.l}" y="${m.t}" width="${plotW}" height="${plotH}"/></clipPath></defs>\n`;
  for (const t of ticks) {
    s += `<line x1="${m.l}" x2="${m.l + plotW}" y1="${f(Y(t))}" y2="${f(Y(t))}" stroke="${GRID}"/>`;
    s += `<text x="${m.l - 8}" y="${f(Y(t) + 4)}" font-size="12" text-anchor="end" fill="${INK2}">${fmtTick(t)}</text>\n`;
  }
  if (50 > y0 && 50 < y1) s += `<line x1="${m.l}" x2="${m.l + plotW}" y1="${f(Y(50))}" y2="${f(Y(50))}" stroke="${INK2}" stroke-dasharray="4 3"/>\n`;
  for (let d = Math.floor(lx0); d <= Math.ceil(lx1); d++) {
    for (let k = 1; k <= 9; k++) {
      const lc = d + Math.log10(k);
      if (lc < lx0 || lc > lx1) continue;
      const x = X(lc);
      s += `<line x1="${f(x)}" x2="${f(x)}" y1="${m.t + plotH}" y2="${m.t + plotH + (k === 1 ? 6 : 3)}" stroke="${INK}"/>`;
      if (k === 1) s += `<text x="${f(x)}" y="${m.t + plotH + 20}" font-size="12" text-anchor="middle" fill="${INK2}">${fmtTick(Math.pow(10, d))}</text>`;
    }
  }
  s += `\n<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + plotH}" stroke="${INK}"/><line x1="${m.l}" x2="${m.l + plotW}" y1="${m.t + plotH}" y2="${m.t + plotH}" stroke="${INK}"/>\n`;
  s += `<text x="${m.l + plotW / 2}" y="${H - 14}" font-size="13" text-anchor="middle" fill="${INK}">Concentration (${esc(o.unit)}, log scale)</text>\n`;
  s += `<text transform="translate(18 ${m.t + plotH / 2}) rotate(-90)" font-size="13" text-anchor="middle" fill="${INK}">Viability (% of untreated control)</text>\n`;

  series.forEach((ser, i) => {
    const col = colorFor(i);
    s += `<g clip-path="url(#plotclip)">`;
    if (ser.fit) {
      const pts: string[] = [];
      for (let k = 0; k <= 200; k++) {
        const lc = lx0 + ((lx1 - lx0) * k) / 200;
        pts.push(`${f(X(lc))},${f(Y(logistic4(lc, ser.fit.params)))}`);
      }
      s += `<polyline points="${pts.join(' ')}" fill="none" stroke="${col}" stroke-width="2"/>`;
    }
    for (const p of ser.points) {
      const x = X(Math.log10(p.conc));
      const tip = `${ser.name}, ${p.conc} ${o.unit}: mean ${p.mean.toFixed(1)}%, SD ${Number.isFinite(p.sd) ? p.sd.toFixed(1) : 'n/a'}, n = ${p.values.length}`;
      if (o.individual) {
        for (const v of p.values) s += `<circle cx="${f(x)}" cy="${f(Y(v))}" r="4" fill="${col}" stroke="#ffffff" stroke-width="1.5"><title>${esc(tip)}</title></circle>`;
      } else {
        if (Number.isFinite(p.sd) && p.sd > 0)
          s += `<line x1="${f(x)}" x2="${f(x)}" y1="${f(Y(p.mean - p.sd))}" y2="${f(Y(p.mean + p.sd))}" stroke="${col}" stroke-width="1.5"/><line x1="${f(x - 4)}" x2="${f(x + 4)}" y1="${f(Y(p.mean - p.sd))}" y2="${f(Y(p.mean - p.sd))}" stroke="${col}" stroke-width="1.5"/><line x1="${f(x - 4)}" x2="${f(x + 4)}" y1="${f(Y(p.mean + p.sd))}" y2="${f(Y(p.mean + p.sd))}" stroke="${col}" stroke-width="1.5"/>`;
        s += `<circle cx="${f(x)}" cy="${f(Y(p.mean))}" r="4.5" fill="${col}" stroke="#ffffff" stroke-width="1.5"><title>${esc(tip)}</title></circle>`;
      }
    }
    s += `</g>\n`;
    const ly = m.t + 6 + i * 40;
    const lx = m.l + plotW + m.r + 10;
    const ic = ser.fit ? `IC50 ${fmtConcPlot(ser.fit.ic50)} ${o.unit}` : 'no fit';
    s += `<line x1="${lx}" x2="${lx + 20}" y1="${ly + 7}" y2="${ly + 7}" stroke="${col}" stroke-width="2"/><circle cx="${lx + 10}" cy="${ly + 7}" r="4" fill="${col}" stroke="#fff" stroke-width="1.5"/>`;
    s += `<text x="${lx + 28}" y="${ly + 11}" font-size="12" fill="${INK}">${esc(ser.name)}</text><text x="${lx + 28}" y="${ly + 26}" font-size="11" fill="${INK2}">${esc(ic)}</text>\n`;
  });
  const note = o.individual ? 'Dots: individual wells' : 'Dots: mean ± SD of replicate wells';
  s += `<text x="${m.l + plotW + m.r + 10}" y="${m.t + 6 + series.length * 40 + 10}" font-size="11" fill="${INK2}">${note}</text>`;
  return s + '</svg>\n';
}

function fmtConcPlot(c: number): string {
  if (!Number.isFinite(c)) return 'n/a';
  if (c >= 1000 || c < 0.01) return c.toExponential(2);
  return Number(c.toPrecision(3)).toString();
}
