// Descriptive statistics and hypothesis tests used by both calculators.
// Distribution functions come from jStat; everything else is implemented here
// and verified against SciPy in tests/stats.test.ts.
import jStat from 'jstat';

export function mean(x: number[]): number {
  if (x.length === 0) return NaN;
  let s = 0;
  for (const v of x) s += v;
  return s / x.length;
}

/** Sample standard deviation (n - 1 denominator). NaN when n < 2. */
export function sd(x: number[]): number {
  const n = x.length;
  if (n < 2) return NaN;
  const m = mean(x);
  let ss = 0;
  for (const v of x) ss += (v - m) ** 2;
  return Math.sqrt(ss / (n - 1));
}

export function sem(x: number[]): number {
  return x.length < 2 ? NaN : sd(x) / Math.sqrt(x.length);
}

export interface Summary {
  n: number;
  mean: number;
  sd: number;
  sem: number;
}

export function summarize(x: number[]): Summary {
  return { n: x.length, mean: mean(x), sd: sd(x), sem: sem(x) };
}

const LANCZOS = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
];

/** log Gamma(x) for x > 0 (Lanczos approximation, g = 7; relative error about 1e-15). */
export function lgamma(x: number): number {
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
  const z = x - 1;
  let a = LANCZOS[0];
  const t = z + 7.5;
  for (let i = 1; i < 9; i++) a += LANCZOS[i] / (z + i);
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(a);
}

/** Continued fraction for the incomplete beta function (modified Lentz method). */
function betacf(x: number, a: number, b: number): number {
  const tiny = 1e-300;
  let c = 1;
  let d = 1 - ((a + b) * x) / (a + 1);
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 1000; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((a + m2 - 1) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (a + b + m) * x) / ((a + m2) * (a + m2 + 1));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return h;
}

/** Regularized incomplete beta function I_x(a, b). */
export function ibeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log1p(-x);
  if (x < (a + 1) / (a + b + 2)) return (Math.exp(lbt) * betacf(x, a, b)) / a;
  return 1 - (Math.exp(lbt) * betacf(1 - x, b, a)) / b;
}

/** Two-sided p-value of a Student t statistic, computed via the incomplete beta function for tail precision. */
export function tTwoSidedP(t: number, df: number): number {
  if (!Number.isFinite(t) || !(df > 0)) return NaN;
  return ibeta(df / (df + t * t), df / 2, 0.5);
}

/** Upper-tail p-value of an F statistic. */
export function fUpperP(f: number, df1: number, df2: number): number {
  if (!Number.isFinite(f) || f < 0 || !(df1 > 0) || !(df2 > 0)) return NaN;
  return ibeta(df2 / (df2 + df1 * f), df2 / 2, df1 / 2);
}

/** Upper quantile of Student's t distribution (p > 0.5): jStat start value refined by Newton steps on our CDF. */
export function tQuantile(p: number, df: number): number {
  let t = jStat.studentt.inv(p, df);
  const logC = lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * Math.log(df * Math.PI);
  for (let i = 0; i < 5; i++) {
    const cdf = t >= 0 ? 1 - tTwoSidedP(t, df) / 2 : tTwoSidedP(t, df) / 2;
    const pdf = Math.exp(logC - ((df + 1) / 2) * Math.log1p((t * t) / df));
    const step = (cdf - p) / pdf;
    t -= step;
    if (Math.abs(step) < 1e-14 * Math.max(1, Math.abs(t))) break;
  }
  return t;
}

export interface WelchResult {
  meanDiff: number; // mean(a) - mean(b)
  t: number;
  df: number;
  p: number;
}

/** Welch's unequal-variance two-sample t-test (two-sided). */
export function welchTTest(a: number[], b: number[]): WelchResult {
  const na = a.length;
  const nb = b.length;
  const va = sd(a) ** 2;
  const vb = sd(b) ** 2;
  const meanDiff = mean(a) - mean(b);
  if (na < 2 || nb < 2) return { meanDiff, t: NaN, df: NaN, p: NaN };
  const sa = va / na;
  const sb = vb / nb;
  const se = Math.sqrt(sa + sb);
  const t = meanDiff / se;
  const df = (sa + sb) ** 2 / (sa ** 2 / (na - 1) + sb ** 2 / (nb - 1));
  return { meanDiff, t, df, p: tTwoSidedP(t, df) };
}

export interface AnovaResult {
  f: number;
  dfBetween: number;
  dfWithin: number;
  p: number;
  mse: number;
}

/** One-way ANOVA across k groups. */
export function oneWayAnova(groups: number[][]): AnovaResult {
  const k = groups.length;
  const all = groups.flat();
  const N = all.length;
  const grand = mean(all);
  let ssb = 0;
  let ssw = 0;
  for (const g of groups) {
    const m = mean(g);
    ssb += g.length * (m - grand) ** 2;
    for (const v of g) ssw += (v - m) ** 2;
  }
  const dfBetween = k - 1;
  const dfWithin = N - k;
  const mse = ssw / dfWithin;
  const f = ssb / dfBetween / mse;
  return { f, dfBetween, dfWithin, p: fUpperP(f, dfBetween, dfWithin), mse };
}

export interface PairwiseResult {
  i: number; // index of the first group
  j: number; // index of the second group
  diff: number; // mean(i) - mean(j)
  stat: number; // q for Tukey, t for Dunnett
  p: number;
}

/** Tukey HSD (Tukey-Kramer for unequal n) for all pairwise comparisons. */
export function tukeyHSD(groups: number[][]): PairwiseResult[] {
  const k = groups.length;
  const { mse, dfWithin } = oneWayAnova(groups);
  const means = groups.map(mean);
  const out: PairwiseResult[] = [];
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const diff = means[i] - means[j];
      const se = Math.sqrt((mse / 2) * (1 / groups[i].length + 1 / groups[j].length));
      const q = Math.abs(diff) / se;
      const p = Math.min(1, Math.max(0, 1 - jStat.tukey.cdf(q, k, dfWithin)));
      out.push({ i, j, diff, stat: q, p });
    }
  }
  return out;
}

function normCdf(x: number): number {
  return jStat.normal.cdf(x, 0, 1);
}

/**
 * P(max_j |T_j| <= c) for Dunnett's many-to-one comparisons.
 * T_j = Z_j / S where Z is multivariate normal with correlation lambda_i * lambda_j
 * and S = sqrt(chi2_df / df). The product correlation structure reduces the
 * multivariate t probability to a two-dimensional integral (Dunnett 1955),
 * evaluated here with composite Simpson rules.
 */
export function dunnettCdf(c: number, lambdas: number[], df: number): number {
  if (c <= 0) return 0;
  const nz = 240;
  const zMax = 8.5;
  const hz = (2 * zMax) / nz;
  const zW: number[] = [];
  const zV: number[] = [];
  for (let a = 0; a <= nz; a++) {
    const z = -zMax + a * hz;
    const w = a === 0 || a === nz ? 1 : a % 2 === 1 ? 4 : 2;
    zV.push(z);
    zW.push((w * hz) / 3 * Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI));
  }
  const scale = lambdas.map((l) => Math.sqrt(1 - l * l));

  const inner = (s: number): number => {
    const cs = c * s;
    let tot = 0;
    for (let a = 0; a <= nz; a++) {
      const z = zV[a];
      let prod = 1;
      for (let m = 0; m < lambdas.length; m++) {
        const lz = lambdas[m] * z;
        prod *= normCdf((cs - lz) / scale[m]) - normCdf((-cs - lz) / scale[m]);
        if (prod < 1e-300) break;
      }
      tot += zW[a] * prod;
    }
    return tot;
  };

  if (!Number.isFinite(df) || df > 1e6) return inner(1);

  // Density of S = sqrt(chi2_df / df).
  const logC = Math.log(2) + (df / 2) * Math.log(df / 2) - lgamma(df / 2);
  const dens = (s: number) => (s <= 0 ? 0 : Math.exp(logC + (df - 1) * Math.log(s) - (df * s * s) / 2));
  const mS = Math.exp(0.5 * Math.log(2 / df) + lgamma((df + 1) / 2) - lgamma(df / 2));
  const sdS = Math.sqrt(Math.max(1 - mS * mS, 1e-12));
  const lo = Math.max(0, mS - 12 * sdS);
  const hi = mS + 14 * sdS;
  const ns = 400;
  const hs = (hi - lo) / ns;
  let total = 0;
  for (let b = 0; b <= ns; b++) {
    const s = lo + b * hs;
    const w = b === 0 || b === ns ? 1 : b % 2 === 1 ? 4 : 2;
    const d = dens(s);
    if (d === 0) continue;
    total += w * d * inner(s);
  }
  return Math.min(1, Math.max(0, (total * hs) / 3));
}

/**
 * Dunnett's test comparing every treatment group with one control group (two-sided, single step).
 * `control` is the index of the control group in `groups`.
 */
export function dunnettTest(groups: number[][], control: number): PairwiseResult[] {
  const { mse, dfWithin } = oneWayAnova(groups);
  const n0 = groups[control].length;
  const m0 = mean(groups[control]);
  const others = groups.map((_, idx) => idx).filter((idx) => idx !== control);
  const lambdas = others.map((idx) => Math.sqrt(groups[idx].length / (groups[idx].length + n0)));
  return others.map((idx) => {
    const g = groups[idx];
    const diff = mean(g) - m0;
    const t = diff / Math.sqrt(mse * (1 / g.length + 1 / n0));
    const p = 1 - dunnettCdf(Math.abs(t), lambdas, dfWithin);
    return { i: idx, j: control, diff, stat: t, p: Math.min(1, Math.max(0, p)) };
  });
}

/** Format a p-value for display. */
export function formatP(p: number): string {
  if (!Number.isFinite(p)) return 'n/a';
  if (p < 0.0001) return '< 0.0001';
  return p.toFixed(4);
}

export function stars(p: number): string {
  if (!Number.isFinite(p)) return '';
  if (p < 0.0001) return '****';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}
