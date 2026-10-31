import { describe, expect, it } from 'vitest';
import { parsePlateMatrix } from '../src/lib/parse';
import { analyzeViability, layoutFromJson, layoutToJson } from '../src/lib/viability';
import { check } from './fit4pl.test';
import { close, readText, ref } from './helpers';

const v = ref.viability;
const { layout } = layoutFromJson(JSON.parse(readText(v.layout)));
const abs = parsePlateMatrix(readText(v.abs));
const refAbs = parsePlateMatrix(readText(v.ref));

describe('viability vs NumPy/SciPy reference', () => {
  it('parses the plate matrices', () => {
    expect(abs.length).toBe(8);
    expect(abs[0].length).toBe(12);
    expect(Number.isNaN(abs[0][11])).toBe(true);
  });

  it('layout JSON round trip', () => {
    const again = layoutFromJson(JSON.parse(JSON.stringify(layoutToJson(layout, 'µM'))));
    expect(again.layout).toEqual(layout);
  });

  it('blank subtraction, reference wavelength and % viability', () => {
    const res = analyzeViability(layout, abs, refAbs);
    close(res.blankMean, v.blankMean, 1e-12, 1e-14, 'blank');
    close(res.controlMean, v.controlMean, 1e-12, 1e-14, 'control');
    for (const [pos, val] of Object.entries(v.viability)) {
      const w = res.wells.find((x) => x.pos === pos)!;
      close(w.viability, val as number, 1e-10, 1e-10, pos);
    }
  });

  it('free 4PL fits per compound', () => {
    const res = analyzeViability(layout, abs, refAbs);
    for (const name of ['Cmpd-A', 'Cmpd-B']) {
      const c = res.compounds.find((x) => x.name === name)!;
      check(c.fit, v.fits[name], name);
      expect(c.flags.some((f) => f.includes('outside the tested range'))).toBe(false);
    }
    // Cmpd-C: only the top of the curve was tested. The free fit is not identifiable and must be flagged.
    const c = res.compounds.find((x) => x.name === 'Cmpd-C')!;
    expect(c.flags.length).toBeGreaterThan(0);
  });

  it('constrained fits (top fixed at 100, bottom fixed at 0)', () => {
    const top = analyzeViability(layout, abs, refAbs, { top: 100 });
    for (const name of ['Cmpd-A', 'Cmpd-B']) check(top.compounds.find((x) => x.name === name)!.fit, v.fits[`${name} (top fixed at 100)`], `${name} top=100`);
    const bot = analyzeViability(layout, abs, refAbs, { bottom: 0 });
    const c = bot.compounds.find((x) => x.name === 'Cmpd-C')!;
    check(c.fit, v.fits['Cmpd-C (bottom fixed at 0)'], 'Cmpd-C bottom=0');
    expect(c.flags.some((f) => f.includes('outside the tested range'))).toBe(true);
  });
});
