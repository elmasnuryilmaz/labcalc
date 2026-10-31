// Four-parameter logistic (4PL) dose-response fit by Levenberg-Marquardt.
//
// Model (x = log10 concentration):
//   y = bottom + (top - bottom) / (1 + 10^(hill * (x - logIC50)))
// The same curve is described by (bottom, top, hill) and (top, bottom, -hill), so
// results are reported in a canonical form with hill > 0: "top" is then the plateau
// at low concentration and "bottom" the plateau at high concentration (an
// increasing curve has bottom > top). The model is identical to LL.4 in the R
// package drc with b = hill, c = bottom, d = top, e = IC50.
import { tQuantile } from './stats';

export interface FitParams {
  bottom: number;
  top: number;
  logIC50: number;
  hill: number;
}

export interface FitResult {
  params: FitParams;
  se: FitParams; // standard errors (NaN when n <= 4)
  ic50: number;
  ic50CI: [number, number]; // 95% CI, back-transformed from log10 scale
  logIC50CI: [number, number];
  absIC50: number; // concentration where the curve crosses 50 (NaN if it never does)
  r2: number;
  ssr: number;
  n: number;
  df: number;
  converged: boolean;
  iterations: number;
  fixed: FixedParams;
}

/** Optional constraints: fix the bottom and/or top plateau at a given value. */
export interface FixedParams {
  bottom?: number;
  top?: number;
}

const LN10 = Math.LN10;

export function logistic4(x: number, p: FitParams): number {
  const e = Math.min(300, Math.max(-300, p.hill * (x - p.logIC50)));
  return p.bottom + (p.top - p.bottom) / (1 + Math.pow(10, e));
}

function fromArr(a: number[]): FitParams {
  return { bottom: a[0], top: a[1], logIC50: a[2], hill: a[3] };
}

/** Jacobian row of the model with respect to (bottom, top, logIC50, hill). */
function jacRow(x: number, a: number[]): number[] {
  const [b, t, c, h] = a;
  const e = Math.min(300, Math.max(-300, h * (x - c)));
  const u = Math.pow(10, e);
  const D = 1 + u;
  const D2 = D * D;
  return [u / D, 1 / D, ((t - b) * u * LN10 * h) / D2, (-(t - b) * u * LN10 * (x - c)) / D2];
}

/** Solve A x = b for a small dense system by Gaussian elimination with partial pivoting. */
export function solve(A: number[][], b: number[]): number[] | null {
  const n = b.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let piv = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
    if (Math.abs(M[piv][col]) < 1e-300) return null;
    [M[col], M[piv]] = [M[piv], M[col]];
    for (let r = col + 1; r < n; r++) {
      const f = M[r][col] / M[col][col];
      for (let k = col; k <= n; k++) M[r][k] -= f * M[col][k];
    }
  }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
    x[r] = s / M[r][r];
  }
  return x;
}

function invert(A: number[][]): number[][] | null {
  const n = A.length;
  const cols: number[][] = [];
  for (let i = 0; i < n; i++) {
    const e = new Array(n).fill(0);
    e[i] = 1;
    const c = solve(A, e);
    if (!c) return null;
    cols.push(c);
  }
  return A.map((_, i) => cols.map((c) => c[i]));
}

function ssrOf(x: number[], y: number[], a: number[]): number {
  const p = fromArr(a);
  let s = 0;
  for (let i = 0; i < x.length; i++) s += (y[i] - logistic4(x[i], p)) ** 2;
  return s;
}

/** Normal equations restricted to the free parameters (indices in `free`). */
function normalMatrix(x: number[], y: number[], a: number[], free: number[]) {
  const k = free.length;
  const JtJ = free.map(() => new Array(k).fill(0));
  const Jtr = new Array(k).fill(0);
  const p = fromArr(a);
  for (let i = 0; i < x.length; i++) {
    const Jf = jacRow(x[i], a);
    const J = free.map((f) => Jf[f]);
    const r = y[i] - logistic4(x[i], p);
    for (let m = 0; m < k; m++) {
      Jtr[m] += J[m] * r;
      for (let q = 0; q < k; q++) JtJ[m][q] += J[m] * J[q];
    }
  }
  return { JtJ, Jtr };
}

function lm(x: number[], y: number[], start: number[], free: number[], maxIter = 1000) {
  let a = [...start];
  let ssr = ssrOf(x, y, a);
  let lambda = 1e-3;
  let converged = false;
  let it = 0;
  for (; it < maxIter; it++) {
    const { JtJ, Jtr } = normalMatrix(x, y, a, free);
    let improved = false;
    for (let inner = 0; inner < 30; inner++) {
      const A = JtJ.map((row, i) => row.map((v, k) => (i === k ? v + lambda * Math.max(v, 1e-12) : v)));
      const delta = solve(A, Jtr);
      if (!delta || delta.some((d) => !Number.isFinite(d))) {
        lambda *= 10;
        continue;
      }
      const cand = [...a];
      free.forEach((f, m) => (cand[f] += delta[m]));
      const s = ssrOf(x, y, cand);
      if (Number.isFinite(s) && s <= ssr) {
        const relChange = (ssr - s) / Math.max(ssr, 1e-300);
        const stepSmall = delta.every((d, m) => Math.abs(d) <= 1e-10 * (Math.abs(a[free[m]]) + 1e-10));
        a = cand;
        ssr = s;
        lambda = Math.max(lambda / 10, 1e-15);
        improved = true;
        if (relChange < 1e-14 || stepSmall) converged = true;
        break;
      }
      lambda *= 10;
      if (lambda > 1e16) break;
    }
    if (!improved) {
      // No downhill step could be found: we are at a (local) minimum.
      converged = true;
      break;
    }
    if (converged) break;
  }
  return { a, ssr, converged, iterations: it + 1 };
}

/** Heuristic starting values from replicate means. */
export function startingValues(x: number[], y: number[]): number[][] {
  const xs = [...new Set(x)].sort((p, q) => p - q);
  const means = xs.map((xv) => {
    const v = y.filter((_, i) => x[i] === xv);
    return v.reduce((s, q) => s + q, 0) / v.length;
  });
  const top = Math.max(...means);
  const bottom = Math.min(...means);
  const half = (top + bottom) / 2;
  let c0 = (xs[0] + xs[xs.length - 1]) / 2;
  for (let i = 0; i < xs.length - 1; i++) {
    const [m1, m2] = [means[i], means[i + 1]];
    if ((m1 - half) * (m2 - half) <= 0 && m1 !== m2) {
      c0 = xs[i] + ((half - m1) / (m2 - m1)) * (xs[i + 1] - xs[i]);
      break;
    }
  }
  const sign = means[0] >= means[means.length - 1] ? 1 : -1;
  const starts: number[][] = [];
  for (const h of [1, 0.5, 2, 4]) {
    for (const c of [c0, xs[0], xs[xs.length - 1]]) starts.push([bottom, top, c, sign * h]);
  }
  return starts;
}

/**
 * Fit the 4PL model to concentration/response pairs (concentrations > 0, linear scale).
 * Several starting points are tried and the lowest residual sum of squares is kept.
 */
export function fit4PL(conc: number[], response: number[], fixed: FixedParams = {}): FitResult | null {
  const pts = conc
    .map((c, i) => [c, response[i]] as const)
    .filter(([c, r]) => c > 0 && Number.isFinite(c) && Number.isFinite(r));
  const x = pts.map(([c]) => Math.log10(c));
  const y = pts.map(([, r]) => r);
  const n = x.length;
  const free = [0, 1, 2, 3].filter((i) => !(i === 0 && fixed.bottom !== undefined) && !(i === 1 && fixed.top !== undefined));
  if (n < free.length || new Set(x).size < 3) return null;

  let best: ReturnType<typeof lm> | null = null;
  for (const s of startingValues(x, y)) {
    if (fixed.bottom !== undefined) s[0] = fixed.bottom;
    if (fixed.top !== undefined) s[1] = fixed.top;
    // With a fixed plateau the sign of the slope is no longer arbitrary: try both.
    for (const st of fixed.bottom !== undefined || fixed.top !== undefined ? [s, [s[0], s[1], s[2], -s[3]]] : [s]) {
      const r = lm(x, y, st, free);
      if (!best || r.ssr < best.ssr - 1e-12 * Math.max(1, best.ssr)) best = r;
    }
  }
  if (!best) return null;
  const a = best.a;
  const df = n - free.length;
  const ybar = y.reduce((s, v) => s + v, 0) / n;
  const sst = y.reduce((s, v) => s + (v - ybar) ** 2, 0);
  const r2 = 1 - best.ssr / sst;

  const seArr = [0, 0, NaN, NaN].map((v, i) => (free.includes(i) ? NaN : v));
  if (df > 0) {
    const { JtJ } = normalMatrix(x, y, a, free);
    const inv = invert(JtJ);
    if (inv) {
      const s2 = best.ssr / df;
      free.forEach((f, m) => (seArr[f] = Math.sqrt(Math.max(0, inv[m][m] * s2))));
    }
  }
  // Canonical form: hill > 0.
  if (a[3] < 0) {
    [a[0], a[1]] = [a[1], a[0]];
    [seArr[0], seArr[1]] = [seArr[1], seArr[0]];
    a[3] = -a[3];
  }
  const params = fromArr(a);
  const se = fromArr(seArr);
  let logIC50CI: [number, number] = [NaN, NaN];
  if (df > 0 && Number.isFinite(se.logIC50)) {
    const tq = tQuantile(0.975, df);
    logIC50CI = [params.logIC50 - tq * se.logIC50, params.logIC50 + tq * se.logIC50];
  }

  // Absolute IC50: concentration where the fitted curve equals 50.
  let absIC50 = NaN;
  const lo = Math.min(params.top, params.bottom);
  const hi = Math.max(params.top, params.bottom);
  if (lo < 50 && hi > 50) {
    const ratio = (params.top - params.bottom) / (50 - params.bottom) - 1;
    if (ratio > 0) absIC50 = Math.pow(10, params.logIC50 + Math.log10(ratio) / params.hill);
  }

  return {
    params,
    se,
    ic50: Math.pow(10, params.logIC50),
    ic50CI: [Math.pow(10, logIC50CI[0]), Math.pow(10, logIC50CI[1])],
    logIC50CI,
    absIC50,
    r2,
    ssr: best.ssr,
    n,
    df,
    converged: best.converged,
    iterations: best.iterations,
    fixed,
  };
}
