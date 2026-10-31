import { describe, expect, it } from 'vitest';
import { guessGroup, parseNumber, parseQpcr, parseWellRanges, splitName } from '../src/lib/parse';
import { analyzeQpcr, buildGeneInfo, buildSampleInfo, defaultSettings, type GeneInfo, type QpcrResult } from '../src/lib/qpcr';
import { close, readText, ref } from './helpers';

function loadLc480() {
  const parsed = parseQpcr(readText(ref.qpcr_lc480.file));
  for (const w of parsed.wells) Object.assign(w, splitName(w.name, '_', 'last'));
  return parsed;
}

function compare(res: QpcrResult, e: any, pfaffl = false) {
  expect(res.perSample.length).toBe(e.perSample.length);
  for (const x of e.perSample) {
    const r = res.perSample.find((p) => p.sample === x.sample && p.gene === x.gene)!;
    expect(r, `${x.sample}/${x.gene}`).toBeDefined();
    close(r.ctTarget, x.ctTarget, 1e-12, 1e-12, 'ct');
    close(r.refCt, x.refCt, 1e-12, 1e-12, 'refCt');
    close(r.dCt, x.dCt, 1e-12, 1e-12, 'dCt');
    close(r.ddCt, x.ddCt, 1e-10, 1e-12, 'ddCt');
    close(r.fc, x.fc, 1e-10, 1e-12, 'FC');
    if (pfaffl) close(r.ratio, x.ratio, 1e-10, 1e-12, 'Pfaffl ratio');
  }
  for (const g of e.groups) {
    const r = res.groupResults.find((q) => q.gene === g.gene && q.group === g.group)!;
    expect(r.dCt.n).toBe(g.n);
    for (const k of ['dCt', 'ddCt', 'fc'] as const)
      for (const s of ['mean', 'sd', 'sem'] as const) close(r[k][s], g[k][s], 1e-10, 1e-12, `${g.gene} ${g.group} ${k} ${s}`);
    close(r.fcOfMean, g.fcOfMean, 1e-10, 0, 'FC of mean ddCt');
    if (pfaffl) for (const s of ['mean', 'sd', 'sem'] as const) close(r.ratio![s], g.ratio[s], 1e-10, 1e-12, `ratio ${s}`);
  }
  const testName: Record<string, string> = { welch: 'Welch', anova: 'One-way ANOVA', tukey: 'Tukey', dunnett: 'Dunnett' };
  for (const s of e.stats) {
    const r = res.stats.find((q) => q.gene === s.gene && q.test.startsWith(testName[s.test]) && (s.test === 'anova' || (q.groupA === s.a && q.groupB === s.b)))!;
    expect(r, `${s.gene} ${s.test} ${s.a}`).toBeDefined();
    if (s.diff !== undefined) close(r.diff, s.diff, 1e-10, 1e-12, 'diff');
    if (s.stat !== undefined) close(r.stat, s.stat, 1e-9, 1e-12, `${s.test} statistic`);
    if (s.df !== undefined) close(r.df as number, s.df, 1e-10, 0, 'df');
    if (s.test === 'dunnett') {
      close(r.p, s.p_exact, 1e-5, 2e-7, 'Dunnett p');
      close(r.p, s.p_scipy, 0, 5e-4, 'Dunnett p vs scipy');
    } else close(r.p, s.p, 1e-5, 1e-7, `${s.test} p`);
  }
}

describe('parsers', () => {
  it('reads the LightCycler 480 Abs Quant export', () => {
    const p = loadLc480();
    expect(p.format).toBe('lc480');
    expect(p.wells.length).toBe(ref.qpcr_lc480.n_wells);
    expect(p.wells[0]).toMatchObject({ pos: 'A1', name: 'Ctrl1_GAPDH', sample: 'Ctrl1', gene: 'GAPDH' });
    expect(p.wells.filter((w) => Number.isNaN(w.ct)).length).toBe(4); // 1 failed well + 3 clean NTCs
  });
  it('numbers, names and ranges', () => {
    expect(parseNumber('23,45')).toBe(23.45);
    expect(Number.isNaN(parseNumber('Undetermined'))).toBe(true);
    expect(splitName('GAPDH Ctrl 1', '', 'first')).toEqual({ gene: 'GAPDH', sample: 'Ctrl 1' });
    expect(guessGroup('TrtA3')).toBe('TrtA');
    expect(guessGroup('KN-4')).toBe('KN');
    expect(parseWellRanges('A1-B3')).toEqual(['A1', 'A2', 'A3', 'B1', 'B2', 'B3']);
  });
  it('handles Include = False and UTF-8 BOM', () => {
    const p = parseQpcr('﻿Include\tColor\tPos\tName\tCp\tConcentration\tStandard\tStatus\nFalse\t255\tA1\tS1\t20.1\t\t0\t\n');
    expect(p.wells[0].excluded).toBe(true);
  });
});

describe('ddCt pipeline vs independent NumPy/SciPy implementation', () => {
  const p = loadLc480();
  const samples = buildSampleInfo(p.wells, {}, guessGroup);
  const genes: Record<string, GeneInfo> = buildGeneInfo(p.wells, {});

  it('auto-detects NTCs, groups and reference genes', () => {
    expect(samples['NTC'].ntc).toBe(true);
    expect(samples['TrtB2'].group).toBe('TrtB');
    expect(genes.GAPDH.reference && genes.ACTB.reference && !genes.IL6.reference).toBe(true);
  });

  it('Livak with Tukey HSD', () => {
    const res = analyzeQpcr(p.wells, samples, genes, { ...defaultSettings, controlGroup: 'Ctrl', multiTest: 'tukey' });
    expect(res.groups).toEqual(['Ctrl', 'TrtA', 'TrtB']);
    compare(res, ref.qpcr_lc480.livak_tukey);
  });

  it('Livak with Dunnett', () => {
    const res = analyzeQpcr(p.wells, samples, genes, { ...defaultSettings, controlGroup: 'Ctrl', multiTest: 'dunnett' });
    compare(res, ref.qpcr_lc480.livak_dunnett);
  });

  it('Pfaffl efficiency-corrected ratio', () => {
    const g2: Record<string, GeneInfo> = {};
    for (const [k, v] of Object.entries(genes)) g2[k] = { ...v, efficiency: ref.qpcr_lc480.efficiency[k] };
    const res = analyzeQpcr(p.wells, samples, g2, { ...defaultSettings, controlGroup: 'Ctrl', method: 'pfaffl' });
    compare(res, ref.qpcr_lc480.pfaffl_tukey, true);
  });

  it('QC flags', () => {
    const res = analyzeQpcr(p.wells, samples, genes, { ...defaultSettings, controlGroup: 'Ctrl' });
    const has = (lvl: string, s: string, g: string, txt: string) => res.qc.some((f) => f.level === lvl && f.sample === s && f.gene === g && f.message.includes(txt));
    expect(has('error', 'TrtB2', 'IL6', 'Technical replicate SD')).toBe(true);
    expect(has('warn', 'TrtA3', 'BCL2', 'No Ct value')).toBe(true);
    expect(has('warn', 'NTC', 'GAPDH', 'NTC amplified')).toBe(true);
    expect(res.qc.filter((f) => f.message.includes('NTC')).length).toBe(1);
  });

  it('excluding a well removes its SD flag', () => {
    const wells = p.wells.map((w) => ({ ...w }));
    const bad = wells.filter((w) => w.sample === 'TrtB2' && w.gene === 'IL6');
    const outlier = bad.reduce((a, b) => (b.ct > a.ct ? b : a));
    outlier.excluded = true;
    const res = analyzeQpcr(wells, samples, genes, { ...defaultSettings, controlGroup: 'Ctrl' });
    expect(res.qc.some((f) => f.sample === 'TrtB2' && f.gene === 'IL6' && f.message.includes('SD'))).toBe(false);
  });
});

describe('two-group CSV with two reference genes (Welch)', () => {
  it('matches the reference', () => {
    const p = parseQpcr(readText(ref.qpcr_two_groups.file));
    expect(p.format).toBe('table');
    const samples = buildSampleInfo(p.wells, {}, guessGroup);
    expect(samples['NTC'].ntc).toBe(true);
    const genes = buildGeneInfo(p.wells, {});
    const res = analyzeQpcr(p.wells, samples, genes, { ...defaultSettings, controlGroup: 'Control' });
    expect(res.references).toEqual(['RPLP0', 'HPRT1']);
    compare(res, ref.qpcr_two_groups.livak);
  });
});
