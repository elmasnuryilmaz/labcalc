import { describe, expect, it } from 'vitest';
import { fit4PL, logistic4 } from '../src/lib/fit4pl';
import { close, ref } from './helpers';

function check(fit: ReturnType<typeof fit4PL>, e: any, label: string) {
  expect(fit, label).not.toBeNull();
  const f = fit!;
  close(f.params.bottom, e.bottom, 1e-5, 1e-5, `${label} bottom`);
  close(f.params.top, e.top, 1e-5, 1e-5, `${label} top`);
  close(f.params.logIC50, e.logIC50, 1e-5, 1e-6, `${label} logIC50`);
  close(f.params.hill, e.hill, 1e-5, 1e-6, `${label} hill`);
  close(f.ic50, e.ic50, 1e-5, 0, `${label} IC50`);
  close(f.se.bottom, e.se[0], 1e-4, 1e-6, `${label} SE bottom`);
  close(f.se.top, e.se[1], 1e-4, 1e-6, `${label} SE top`);
  close(f.se.logIC50, e.se[2], 1e-4, 1e-8, `${label} SE logIC50`);
  close(f.se.hill, e.se[3], 1e-4, 1e-8, `${label} SE hill`);
  close(f.ic50CI[0], e.ic50CI[0], 1e-4, 0, `${label} CI low`);
  close(f.ic50CI[1], e.ic50CI[1], 1e-4, 0, `${label} CI high`);
  close(f.r2, e.r2, 1e-8, 1e-10, `${label} R2`);
  close(f.ssr, e.ssr, 1e-7, 1e-10, `${label} SSR`);
  close(f.absIC50, e.absIC50, 1e-5, 0, `${label} absolute IC50`);
  expect(f.df).toBe(e.df);
}

describe('4PL fit vs scipy.optimize.curve_fit', () => {
  it('increasing curve (canonical form, hill > 0, bottom > top)', () => {
    const c = ref.fit_cases[0];
    check(fit4PL(c.conc, c.y), c, 'increasing');
  });
  it('five points, one residual degree of freedom', () => {
    const c = ref.fit_cases[1];
    check(fit4PL(c.conc, c.y), c, 'n=5');
  });
  it('top and bottom fixed (two-parameter fit)', () => {
    const c = ref.fit_cases[2];
    check(fit4PL(c.conc, c.y, { bottom: c.fixBottom, top: c.fixTop }), c, 'fixed');
  });
  it('returns null with too few points', () => {
    expect(fit4PL([1, 2, 3], [90, 50, 10])).toBeNull();
  });
  it('recovers exact parameters from noise-free data', () => {
    const p = { bottom: 3, top: 97, logIC50: -0.5, hill: 1.7 };
    const conc = [0.001, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 100];
    const f = fit4PL(conc, conc.map((c) => logistic4(Math.log10(c), p)))!;
    close(f.params.logIC50, p.logIC50, 1e-8, 1e-9);
    close(f.params.hill, p.hill, 1e-7);
    close(f.r2, 1, 0, 1e-12);
  });
});

export { check };
