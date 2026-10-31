// Relative quantification for SYBR Green qPCR: Livak 2^-ddCt and Pfaffl
// efficiency-corrected ratios, quality-control flags and group statistics.
import type { Well } from './parse';
import {
  dunnettTest,
  mean,
  oneWayAnova,
  sd,
  summarize,
  tukeyHSD,
  welchTTest,
  type Summary,
} from './stats';

export interface SampleInfo {
  group: string;
  ntc: boolean;
}

export interface GeneInfo {
  reference: boolean;
  efficiency: number | null; // amplification efficiency in percent (100 = perfect doubling)
}

export interface QpcrSettings {
  sdThreshold: number; // technical replicate SD above this gives a warning
  sdMax: number; // technical replicate SD above this is flagged as an error
  ctMax: number; // Ct above this is flagged as late amplification
  controlGroup: string;
  method: 'livak' | 'pfaffl';
  multiTest: 'tukey' | 'dunnett';
}

export const defaultSettings: QpcrSettings = {
  sdThreshold: 0.3,
  sdMax: 0.5,
  ctMax: 35,
  controlGroup: '',
  method: 'livak',
  multiTest: 'tukey',
};

export interface QcFlag {
  level: 'warn' | 'error';
  sample: string;
  gene: string;
  wells: string[];
  message: string;
}

export interface TechRep {
  sample: string;
  gene: string;
  wellIds: number[];
  cts: number[];
  ctMean: number;
  ctSd: number;
}

export interface SampleResult {
  sample: string;
  group: string;
  gene: string;
  ctTarget: number;
  ctTargetSd: number;
  refCt: number; // mean of reference-gene Ct (log of their geometric mean quantity)
  dCt: number;
  ddCt: number;
  fc: number; // 2^-ddCt
  dCqE: number; // efficiency-weighted dCt (NaN unless efficiencies are known)
  log2Ratio: number; // Pfaffl log2 ratio
  ratio: number; // Pfaffl ratio
}

export interface GroupResult {
  gene: string;
  group: string;
  dCt: Summary;
  ddCt: Summary;
  fc: Summary;
  fcOfMean: number; // 2^-(mean ddCt)
  ratio: Summary | null;
  log2Ratio: Summary | null;
}

export interface StatRow {
  gene: string;
  test: string;
  comparison: string;
  groupA: string;
  groupB: string;
  diff: number; // difference in mean dCt (A - B); NaN for ANOVA rows
  stat: number;
  df: number | string;
  p: number;
}

export interface QpcrResult {
  techReps: TechRep[];
  qc: QcFlag[];
  targets: string[];
  references: string[];
  groups: string[]; // control group first
  perSample: SampleResult[];
  groupResults: GroupResult[];
  stats: StatRow[];
  notes: string[];
}

const key = (s: string, g: string) => `${s}\u0000${g}`;

/** Amplification factor from percent efficiency: 100 % -> 2. */
export const ampFactor = (effPercent: number) => 1 + effPercent / 100;

export function analyzeQpcr(
  wells: Well[],
  samples: Record<string, SampleInfo>,
  genes: Record<string, GeneInfo>,
  settings: QpcrSettings,
): QpcrResult {
  const qc: QcFlag[] = [];
  const notes: string[] = [];
  const isNtc = (s: string) => samples[s]?.ntc ?? false;

  const unassigned = wells.filter((w) => !w.excluded && (!w.sample || !w.gene));
  if (unassigned.length)
    qc.push({ level: 'warn', sample: '', gene: '', wells: unassigned.map((w) => w.pos || `#${w.id + 1}`), message: `${unassigned.length} well(s) without sample or gene are ignored.` });

  // Well-level flags and technical replicate grouping.
  const reps = new Map<string, TechRep>();
  for (const w of wells) {
    if (w.excluded || !w.sample || !w.gene) continue;
    const label = w.pos || `#${w.id + 1}`;
    if (isNtc(w.sample)) {
      if (Number.isFinite(w.ct))
        qc.push({
          level: w.ct < settings.ctMax ? 'error' : 'warn',
          sample: w.sample,
          gene: w.gene,
          wells: [label],
          message: `NTC amplified (Ct ${w.ct.toFixed(2)})${w.ct < settings.ctMax ? '' : ', probably primer dimers'}. Check the melt curve.`,
        });
      continue;
    }
    if (!Number.isFinite(w.ct)) {
      qc.push({ level: 'warn', sample: w.sample, gene: w.gene, wells: [label], message: w.ctRaw ? `Undetermined Ct ("${w.ctRaw}"); well ignored.` : 'No Ct value (no amplification or missing); well ignored.' });
      continue;
    }
    if (w.ct > settings.ctMax)
      qc.push({ level: 'warn', sample: w.sample, gene: w.gene, wells: [label], message: `Ct ${w.ct.toFixed(2)} above ${settings.ctMax} (late amplification, low reliability).` });
    const k = key(w.sample, w.gene);
    let r = reps.get(k);
    if (!r) {
      r = { sample: w.sample, gene: w.gene, wellIds: [], cts: [], ctMean: NaN, ctSd: NaN };
      reps.set(k, r);
    }
    r.wellIds.push(w.id);
    r.cts.push(w.ct);
  }
  const posOf = (id: number) => {
    const w = wells.find((x) => x.id === id)!;
    return w.pos || `#${w.id + 1}`;
  };
  for (const r of reps.values()) {
    r.ctMean = mean(r.cts);
    r.ctSd = sd(r.cts);
    if (r.cts.length >= 2 && r.ctSd > settings.sdThreshold) {
      const err = r.ctSd > settings.sdMax;
      qc.push({
        level: err ? 'error' : 'warn',
        sample: r.sample,
        gene: r.gene,
        wells: r.wellIds.map(posOf),
        message: `Technical replicate SD ${r.ctSd.toFixed(2)} > ${(err ? settings.sdMax : settings.sdThreshold).toFixed(2)}.${err ? ' Consider excluding the outlier well.' : ''}`,
      });
    }
    if (r.cts.length === 1) qc.push({ level: 'warn', sample: r.sample, gene: r.gene, wells: r.wellIds.map(posOf), message: 'Only one valid technical replicate.' });
  }

  const allGenes = [...new Set(wells.filter((w) => w.gene).map((w) => w.gene))];
  const references = allGenes.filter((g) => genes[g]?.reference);
  const targets = allGenes.filter((g) => !genes[g]?.reference);
  const sampleNames = [...new Set(wells.filter((w) => w.sample && !isNtc(w.sample)).map((w) => w.sample))];
  const groupsSeen: string[] = [];
  for (const s of sampleNames) {
    const g = samples[s]?.group ?? '';
    if (g && !groupsSeen.includes(g)) groupsSeen.push(g);
  }
  const control = groupsSeen.includes(settings.controlGroup) ? settings.controlGroup : groupsSeen[0] ?? '';
  const groups = control ? [control, ...groupsSeen.filter((g) => g !== control)] : [];

  const result: QpcrResult = { techReps: [...reps.values()], qc, targets, references, groups, perSample: [], groupResults: [], stats: [], notes };
  if (references.length === 0) {
    notes.push('Mark at least one gene as a reference gene to calculate dCt.');
    return result;
  }
  if (targets.length === 0) {
    notes.push('All genes are marked as reference genes; there is no target gene.');
    return result;
  }
  const effKnown = (g: string) => genes[g]?.efficiency != null && Number.isFinite(genes[g]!.efficiency!);
  const log2E = (g: string) => Math.log2(ampFactor(genes[g]!.efficiency!));
  const refsHaveEff = references.every(effKnown);

  for (const s of sampleNames) {
    const group = samples[s]?.group ?? '';
    if (!group) continue;
    const refCts = references.map((g) => reps.get(key(s, g))?.ctMean ?? NaN);
    const missingRef = references.filter((_, i) => !Number.isFinite(refCts[i]));
    if (missingRef.length) {
      qc.push({ level: 'error', sample: s, gene: missingRef.join(', '), wells: [], message: `No valid Ct for reference gene(s) ${missingRef.join(', ')}; sample dropped.` });
      continue;
    }
    const refCt = mean(refCts);
    for (const g of targets) {
      const r = reps.get(key(s, g));
      if (!r) {
        qc.push({ level: 'error', sample: s, gene: g, wells: [], message: `No valid Ct for ${g}; sample dropped for this gene.` });
        continue;
      }
      const dCt = r.ctMean - refCt;
      let dCqE = NaN;
      if (effKnown(g) && refsHaveEff) dCqE = log2E(g) * r.ctMean - mean(references.map((rg, i) => log2E(rg) * refCts[i]));
      result.perSample.push({ sample: s, group, gene: g, ctTarget: r.ctMean, ctTargetSd: r.ctSd, refCt, dCt, ddCt: NaN, fc: NaN, dCqE, log2Ratio: NaN, ratio: NaN });
    }
  }

  // Calibrate to the control group.
  for (const g of targets) {
    const rows = result.perSample.filter((r) => r.gene === g);
    const ctrl = rows.filter((r) => r.group === control);
    if (ctrl.length === 0) {
      notes.push(`${g}: no control-group samples, so ddCt cannot be calculated.`);
      continue;
    }
    const ctrlMean = mean(ctrl.map((r) => r.dCt));
    const ctrlMeanE = mean(ctrl.map((r) => r.dCqE));
    for (const r of rows) {
      r.ddCt = r.dCt - ctrlMean;
      r.fc = Math.pow(2, -r.ddCt);
      if (Number.isFinite(r.dCqE)) {
        r.log2Ratio = -(r.dCqE - ctrlMeanE);
        r.ratio = Math.pow(2, r.log2Ratio);
      }
    }
  }
  if (settings.method === 'pfaffl' && !targets.every((g) => effKnown(g) && refsHaveEff))
    notes.push('Pfaffl ratios need efficiencies for the target and every reference gene; genes without them fall back to Livak for plots and statistics.');

  // Group summaries.
  for (const g of targets) {
    for (const grp of groups) {
      const rows = result.perSample.filter((r) => r.gene === g && r.group === grp && Number.isFinite(r.ddCt));
      if (rows.length === 0) continue;
      const hasE = rows.every((r) => Number.isFinite(r.ratio));
      result.groupResults.push({
        gene: g,
        group: grp,
        dCt: summarize(rows.map((r) => r.dCt)),
        ddCt: summarize(rows.map((r) => r.ddCt)),
        fc: summarize(rows.map((r) => r.fc)),
        fcOfMean: Math.pow(2, -mean(rows.map((r) => r.ddCt))),
        ratio: hasE ? summarize(rows.map((r) => r.ratio)) : null,
        log2Ratio: hasE ? summarize(rows.map((r) => r.log2Ratio)) : null,
      });
    }
  }

  // Statistics on dCt (or efficiency-weighted dCt for Pfaffl), never on fold change.
  for (const g of targets) {
    const usePfaffl = settings.method === 'pfaffl' && effKnown(g) && refsHaveEff;
    const val = (r: SampleResult) => (usePfaffl ? r.dCqE : r.dCt);
    const label = usePfaffl ? ' (efficiency-weighted dCt)' : ' (dCt)';
    const present = groups.filter((grp) => result.perSample.some((r) => r.gene === g && r.group === grp && Number.isFinite(r.ddCt)));
    const data = present.map((grp) => result.perSample.filter((r) => r.gene === g && r.group === grp && Number.isFinite(r.ddCt)).map(val));
    if (present.length < 2 || present[0] !== control) continue;
    if (data.some((d) => d.length < 2)) {
      notes.push(`${g}: statistics need at least 2 biological replicates per group.`);
      continue;
    }
    if (present.length === 2) {
      const w = welchTTest(data[1], data[0]);
      result.stats.push({ gene: g, test: 'Welch t-test' + label, comparison: `${present[1]} vs ${present[0]}`, groupA: present[1], groupB: present[0], diff: w.meanDiff, stat: w.t, df: w.df, p: w.p });
    } else {
      const a = oneWayAnova(data);
      result.stats.push({ gene: g, test: 'One-way ANOVA' + label, comparison: 'all groups', groupA: '', groupB: '', diff: NaN, stat: a.f, df: `${a.dfBetween}, ${a.dfWithin}`, p: a.p });
      if (settings.multiTest === 'tukey') {
        for (const pr of tukeyHSD(data))
          result.stats.push({ gene: g, test: 'Tukey HSD', comparison: `${present[pr.j]} vs ${present[pr.i]}`, groupA: present[pr.j], groupB: present[pr.i], diff: -pr.diff, stat: pr.stat, df: a.dfWithin, p: pr.p });
      } else {
        for (const pr of dunnettTest(data, 0))
          result.stats.push({ gene: g, test: 'Dunnett vs control', comparison: `${present[pr.i]} vs ${present[0]}`, groupA: present[pr.i], groupB: present[0], diff: pr.diff, stat: pr.stat, df: a.dfWithin, p: pr.p });
      }
    }
  }
  return result;
}

const NTC_RE = /^(ntc|h2o|water|nrt|neg)/i;

/**
 * Build (or extend) the sample table from the wells: the group comes from a group
 * column when the file has one, otherwise it is guessed from the sample name;
 * samples called NTC, H2O, water, NRT or neg* are treated as no-template controls.
 */
export function buildSampleInfo(wells: Well[], existing: Record<string, SampleInfo>, guess: (s: string) => string): Record<string, SampleInfo> {
  const out: Record<string, SampleInfo> = {};
  for (const w of wells) {
    if (!w.sample || out[w.sample]) continue;
    if (existing[w.sample]) {
      out[w.sample] = existing[w.sample];
      continue;
    }
    const ntc = NTC_RE.test(w.sample) || /^ntc$/i.test(w.group);
    out[w.sample] = { group: ntc ? '' : w.group || guess(w.sample), ntc };
  }
  return out;
}

export function buildGeneInfo(wells: Well[], existing: Record<string, GeneInfo>): Record<string, GeneInfo> {
  const out: Record<string, GeneInfo> = {};
  for (const w of wells) {
    if (!w.gene || out[w.gene]) continue;
    out[w.gene] = existing[w.gene] ?? { reference: COMMON_REFS.test(w.gene), efficiency: null };
  }
  return out;
}

/** Commonly used reference genes, pre-ticked as references (the user can change this). */
const COMMON_REFS = /^(gapdh|actb|b2m|rplp0|hprt1?|tbp|ppia|ywhaz|ubc|sdha|gusb|18s|rn18s|hmbs|tubb|eef1a1|pgk1|rpl13a|rps18)$/i;
