import { describe, expect, it } from 'vitest';
import jStat from 'jstat';
import { dunnettTest, fUpperP, mean, oneWayAnova, sd, tQuantile, tTwoSidedP, tukeyHSD, welchTTest } from '../src/lib/stats';
import { close, ref } from './helpers';

describe('distribution functions vs SciPy', () => {
  it('two-sided t p-values', () => {
    for (const c of ref.dist.t_two_sided) close(tTwoSidedP(c.t, c.df), c.p, 1e-11, 1e-15, `t=${c.t} df=${c.df}`);
  });
  it('F upper-tail p-values', () => {
    for (const c of ref.dist.f_upper) close(fUpperP(c.f, c.df1, c.df2), c.p, 1e-11, 1e-15, `F=${c.f}`);
  });
  it('t quantiles', () => {
    for (const c of ref.dist.t_quantile) close(tQuantile(c.p, c.df), c.q, 1e-11, 0, `df=${c.df}`);
  });
  it('studentized range CDF', () => {
    for (const c of ref.dist.tukey_cdf) close(jStat.tukey.cdf(c.q, c.k, c.df), c.cdf, 1e-7, 1e-9, `q=${c.q} k=${c.k} df=${c.df}`);
  });
});

describe('descriptive statistics', () => {
  it('mean and sample SD', () => {
    expect(mean([1, 2, 3, 4])).toBe(2.5);
    close(sd([1, 2, 3, 4]), Math.sqrt(5 / 3), 1e-14);
    expect(Number.isNaN(sd([1]))).toBe(true);
  });
});

describe('Welch t-test vs scipy.stats.ttest_ind(equal_var=False)', () => {
  it.each(ref.welch.map((c: unknown, i: number) => [i, c]) as [number, any][])('case %i', (_i, c) => {
    const r = welchTTest(c.a, c.b);
    close(r.t, c.t, 1e-10, 0, 't');
    close(r.df, c.df, 1e-10, 0, 'df');
    close(r.p, c.p, 1e-10, 1e-14, 'p');
  });
});

describe('ANOVA, Tukey HSD and Dunnett vs SciPy', () => {
  it.each(ref.multi.map((c: unknown, i: number) => [i, c]) as [number, any][])('case %i', (_i, c) => {
    const a = oneWayAnova(c.groups);
    close(a.f, c.anova.f, 1e-10, 0, 'F');
    expect(a.dfBetween).toBe(c.anova.df1);
    expect(a.dfWithin).toBe(c.anova.df2);
    close(a.p, c.anova.p, 1e-10, 1e-14, 'ANOVA p');

    const tk = tukeyHSD(c.groups);
    expect(tk.length).toBe(c.tukey.length);
    tk.forEach((r, idx) => {
      const e = c.tukey[idx];
      expect([r.i, r.j]).toEqual([e.i, e.j]);
      close(r.diff, e.diff, 1e-10, 1e-12, 'Tukey diff');
      close(r.stat, e.q, 1e-10, 0, 'Tukey q');
      close(r.p, e.p, 1e-5, 1e-7, `Tukey p ${e.i}-${e.j}`);
    });

    const dn = dunnettTest(c.groups, 0);
    dn.forEach((r, idx) => {
      const e = c.dunnett[idx];
      close(r.stat, e.t, 1e-10, 0, 'Dunnett t');
      // The exact value is a nested adaptive quadrature in SciPy; scipy.stats.dunnett itself is Monte Carlo based.
      close(r.p, e.p_exact, 1e-5, 2e-7, 'Dunnett p (quadrature)');
      close(r.p, e.p_scipy, 0, 5e-4, 'Dunnett p (scipy.stats.dunnett)');
    });
  });

  it('Dunnett with a single comparison equals the pooled two-sample t-test', () => {
    const g = [[1.1, 2.3, 1.9, 2.4], [3.2, 2.8, 3.9]];
    const [r] = dunnettTest(g, 0);
    const a = oneWayAnova(g);
    const t = (mean(g[1]) - mean(g[0])) / Math.sqrt(a.mse * (1 / 3 + 1 / 4));
    close(r.stat, t, 1e-12);
    close(r.p, tTwoSidedP(t, a.dfWithin), 1e-6, 1e-8);
  });
});
