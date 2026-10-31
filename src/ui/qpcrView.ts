import lc480Example from '../../examples/qpcr_lc480_synthetic.txt?raw';
import twoGroupExample from '../../examples/qpcr_two_groups_synthetic.csv?raw';
import { download, downloadPng, downloadSvg, toCsv } from '../lib/export';
import { guessGroup, parseQpcr, parseWellRanges, splitName, normPos, type Well } from '../lib/parse';
import { barPlot, type BarSeries } from '../lib/plot';
import {
  analyzeQpcr,
  ampFactor,
  buildGeneInfo,
  buildSampleInfo,
  defaultSettings,
  type GeneInfo,
  type QpcrResult,
  type QpcrSettings,
  type SampleInfo,
} from '../lib/qpcr';
import { formatP } from '../lib/stats';
import { card, fmt, h, readFileText, rerender, table } from './dom';

interface PlotOpts {
  log2: boolean;
  error: 'sd' | 'sem';
  stars: boolean;
  points: boolean;
  title: string;
  hidden: string[]; // target genes not shown
}

interface State {
  source: string;
  format: string;
  parseWarnings: string[];
  error: string;
  wells: Well[];
  samples: Record<string, SampleInfo>;
  genes: Record<string, GeneInfo>;
  settings: QpcrSettings;
  plot: PlotOpts;
  split: { delim: string; genePart: 'first' | 'last' };
  bulk: { wells: string; field: 'sample' | 'gene'; value: string };
}

const st: State = {
  source: '',
  format: '',
  parseWarnings: [],
  error: '',
  wells: [],
  samples: {},
  genes: {},
  settings: { ...defaultSettings },
  plot: { log2: false, error: 'sd', stars: true, points: true, title: '', hidden: [] },
  split: { delim: '_', genePart: 'last' },
  bulk: { wells: '', field: 'gene', value: '' },
};

let root: HTMLElement;
let pasteText = '';

function refreshTables() {
  st.samples = buildSampleInfo(st.wells, st.samples, guessGroup);
  st.genes = buildGeneInfo(st.wells, st.genes);
  const groups = [...new Set(Object.values(st.samples).filter((s) => !s.ntc && s.group).map((s) => s.group))];
  if (!groups.includes(st.settings.controlGroup))
    st.settings.controlGroup = groups.find((g) => /^(ctrl|control|untreated|vehicle|dmso|mock|con|wt|sham)/i.test(g)) ?? groups[0] ?? '';
}

function load(text: string, source: string) {
  try {
    const p = parseQpcr(text);
    st.wells = p.wells;
    st.format = p.format === 'lc480' ? 'LightCycler 480 export' : 'table (CSV/TSV)';
    st.parseWarnings = p.warnings;
    st.error = '';
    st.source = source;
    st.samples = {};
    st.genes = {};
    st.settings.controlGroup = '';
    // LightCycler names usually hold "sample_gene"; split them right away when they do.
    if (st.wells.some((w) => !w.gene) && st.wells.filter((w) => w.name.includes('_')).length > st.wells.length / 2)
      applySplit();
    else for (const w of st.wells) if (!w.sample) w.sample = w.name;
    refreshTables();
  } catch (e) {
    st.error = (e as Error).message;
    st.wells = [];
  }
  render();
}

function applySplit() {
  for (const w of st.wells) {
    const r = splitName(w.name, st.split.delim, st.split.genePart);
    w.sample = r.sample;
    w.gene = r.gene;
  }
}

function update(fn: () => void, tables = false) {
  fn();
  if (tables) refreshTables();
  // Defer so that focus has moved to the next field before re-rendering.
  setTimeout(render, 0);
}

export function mountQpcr(el: HTMLElement) {
  root = el;
  render();
}

function render() {
  rerender(root, build);
}

function build(): Node {
  const res = st.wells.length ? analyzeQpcr(st.wells, st.samples, st.genes, st.settings) : null;
  return h('div', null, loadCard(), st.wells.length ? [wellCard(res!), mappingCard(), settingsCard(res!), qcCard(res!), resultsCard(res!), figureCard(res!)] : null);
}

function loadCard() {
  const file = h('input', {
    type: 'file',
    accept: '.txt,.csv,.tsv,.tab',
    onchange: async (e: Event) => {
      const f = (e.target as HTMLInputElement).files?.[0];
      if (f) load(await readFileText(f), f.name);
    },
  });
  return card(
    '1. Load qPCR data',
    h('p', { class: 'hint' }, 'Upload or paste a Roche LightCycler 480 text export (Abs Quant / 2nd Derivative Max: Include, Color, Pos, Name, Cp, ...) or any CSV/TSV table with a header row and a Ct, Cp or Cq column. Optional columns: Well, Sample, Gene/Target, Group. SYBR Green data; one row per well.'),
    h('div', { class: 'row' }, h('label', { class: 'btn' }, 'Choose file', file), h('button', { onclick: () => load(lc480Example, 'Synthetic LightCycler 480 example') }, 'Load synthetic LC480 example (3 groups)'), h('button', { onclick: () => load(twoGroupExample, 'Synthetic two-group CSV example') }, 'Load synthetic CSV example (2 groups)')),
    h('details', null, h('summary', null, 'Paste data instead'), h('textarea', { rows: 6, placeholder: 'Paste the export here', oninput: (e: Event) => (pasteText = (e.target as HTMLTextAreaElement).value) }), h('button', { onclick: () => load(pasteText, 'pasted text') }, 'Read pasted data')),
    st.error ? h('p', { class: 'flag error' }, st.error) : null,
    st.wells.length ? h('p', { class: 'ok' }, `Loaded ${st.wells.length} wells from ${st.source} (${st.format}).${st.source.startsWith('Synthetic') ? ' These are simulated data for demonstration only.' : ''}`) : null,
    st.parseWarnings.map((w) => h('p', { class: 'flag warn' }, w)),
  );
}

function wellCard(res: QpcrResult) {
  const flagged = new Map<string, string>();
  for (const f of res.qc) for (const w of f.wells) if (!flagged.has(w) || f.level === 'error') flagged.set(w, f.level);
  const rows = st.wells.map((w) => {
    const label = w.pos || `#${w.id + 1}`;
    const lvl = flagged.get(label);
    return [
      label,
      w.name,
      Number.isFinite(w.ct) ? w.ct.toFixed(2) : h('span', { class: 'muted' }, w.ctRaw || 'missing'),
      h('input', { value: w.sample, 'data-k': `s${w.id}`, onchange: (e: Event) => update(() => (w.sample = (e.target as HTMLInputElement).value.trim()), true) }),
      h('input', { value: w.gene, 'data-k': `g${w.id}`, onchange: (e: Event) => update(() => (w.gene = (e.target as HTMLInputElement).value.trim()), true) }),
      h('input', { type: 'checkbox', checked: w.excluded, 'data-k': `x${w.id}`, onchange: (e: Event) => update(() => (w.excluded = (e.target as HTMLInputElement).checked)) }),
      lvl ? h('span', { class: `badge ${lvl}` }, lvl === 'error' ? 'check' : 'note') : '',
    ];
  });
  return card(
    '2. Assign wells to samples and genes',
    h('div', { class: 'row' },
      h('span', null, 'Split the Name column at'),
      h('input', { class: 'short', value: st.split.delim, title: 'Leave empty to split at spaces', onchange: (e: Event) => (st.split.delim = (e.target as HTMLInputElement).value) }),
      h('span', null, 'gene name is the'),
      h('select', { onchange: (e: Event) => (st.split.genePart = (e.target as HTMLSelectElement).value as 'first' | 'last') }, h('option', { value: 'last', selected: st.split.genePart === 'last' }, 'last part'), h('option', { value: 'first', selected: st.split.genePart === 'first' }, 'first part')),
      h('button', { onclick: () => update(applySplit, true) }, 'Split names'),
      h('button', { onclick: () => update(() => st.wells.forEach((w) => ((w.sample = w.name), (w.gene = ''))), true) }, 'Use names as samples'),
    ),
    h('div', { class: 'row' },
      h('span', null, 'Set'),
      h('select', { onchange: (e: Event) => (st.bulk.field = (e.target as HTMLSelectElement).value as 'sample' | 'gene') }, h('option', { value: 'gene', selected: st.bulk.field === 'gene' }, 'gene'), h('option', { value: 'sample', selected: st.bulk.field === 'sample' }, 'sample')),
      h('span', null, 'of wells'),
      h('input', { placeholder: 'e.g. A1-A12, B1:B6', value: st.bulk.wells, onchange: (e: Event) => (st.bulk.wells = (e.target as HTMLInputElement).value) }),
      h('span', null, 'to'),
      h('input', { placeholder: 'value', value: st.bulk.value, onchange: (e: Event) => (st.bulk.value = (e.target as HTMLInputElement).value.trim()) }),
      h('button', {
        onclick: () => {
          try {
            const set = new Set(parseWellRanges(st.bulk.wells));
            update(() => st.wells.filter((w) => set.has(normPos(w.pos))).forEach((w) => (w[st.bulk.field] = st.bulk.value)), true);
          } catch (e) {
            alert((e as Error).message);
          }
        },
      }, 'Apply'),
    ),
    h('div', { class: 'scroll', 'data-scroll': 'wells' }, table(['Well', 'Name', 'Ct', 'Sample', 'Gene', 'Exclude', 'QC'], rows, 'wells')),
  );
}

function mappingCard() {
  const groups = [...new Set(Object.values(st.samples).map((s) => s.group).filter(Boolean))];
  const sampleRows = Object.entries(st.samples).map(([name, s]) => [
    name,
    h('input', { value: s.group, list: 'grouplist', disabled: s.ntc, 'data-k': `grp-${name}`, onchange: (e: Event) => update(() => (s.group = (e.target as HTMLInputElement).value.trim()), true) }),
    h('input', { type: 'checkbox', checked: s.ntc, 'data-k': `ntc-${name}`, onchange: (e: Event) => update(() => (s.ntc = (e.target as HTMLInputElement).checked), true) }),
  ]);
  const geneRows = Object.entries(st.genes).map(([name, g]) => [
    h('i', null, name),
    h('input', { type: 'checkbox', checked: g.reference, 'data-k': `ref-${name}`, onchange: (e: Event) => update(() => (g.reference = (e.target as HTMLInputElement).checked)) }),
    h('input', {
      class: 'short',
      type: 'number',
      step: '0.1',
      min: '50',
      max: '150',
      placeholder: '%',
      value: g.efficiency ?? '',
      'data-k': `eff-${name}`,
      onchange: (e: Event) =>
        update(() => {
          const v = parseFloat((e.target as HTMLInputElement).value);
          g.efficiency = Number.isFinite(v) ? v : null;
        }),
    }),
    g.efficiency != null ? `E = ${ampFactor(g.efficiency).toFixed(3)}${g.efficiency < 90 || g.efficiency > 110 ? ' (outside 90 to 110 %)' : ''}` : '',
  ]);
  return card(
    '3. Samples, groups and genes',
    h('datalist', { id: 'grouplist' }, groups.map((g) => h('option', { value: g }))),
    h('div', { class: 'cols' },
      h('div', null, h('h3', null, 'Biological samples'), h('p', { class: 'hint' }, 'Technical replicates (wells with the same sample and gene) are averaged first. Tick NTC for no-template controls.'), h('div', { class: 'scroll short', 'data-scroll': 'samples' }, table(['Sample', 'Group', 'NTC'], sampleRows))),
      h('div', null, h('h3', null, 'Genes'), h('p', { class: 'hint' }, 'Tick the reference gene(s); with several, their Ct values are averaged (geometric mean of quantities). Efficiencies (%) are optional and enable the Pfaffl ratio.'), table(['Gene', 'Reference', 'Efficiency %', ''], geneRows)),
    ),
  );
}

function settingsCard(res: QpcrResult) {
  const s = st.settings;
  const num = (k: 'sdThreshold' | 'sdMax' | 'ctMax', step: string) =>
    h('input', { type: 'number', class: 'short', step, value: s[k], 'data-k': `set-${k}`, onchange: (e: Event) => update(() => (s[k] = parseFloat((e.target as HTMLInputElement).value) || defaultSettings[k])) });
  return card(
    '4. Analysis settings',
    h('div', { class: 'grid2' },
      h('label', null, 'Control (calibrator) group'),
      h('select', { onchange: (e: Event) => update(() => (s.controlGroup = (e.target as HTMLSelectElement).value)) }, res.groups.map((g) => h('option', { value: g, selected: g === s.controlGroup }, g))),
      h('label', null, 'Method for figure and statistics'),
      h('select', { onchange: (e: Event) => update(() => (s.method = (e.target as HTMLSelectElement).value as 'livak' | 'pfaffl')) }, h('option', { value: 'livak', selected: s.method === 'livak' }, 'Livak 2^-ΔΔCt (assumes 100 % efficiency)'), h('option', { value: 'pfaffl', selected: s.method === 'pfaffl' }, 'Pfaffl (efficiency-corrected)')),
      h('label', null, 'More than two groups'),
      h('select', { onchange: (e: Event) => update(() => (s.multiTest = (e.target as HTMLSelectElement).value as 'tukey' | 'dunnett')) }, h('option', { value: 'tukey', selected: s.multiTest === 'tukey' }, 'One-way ANOVA + Tukey HSD (all pairs)'), h('option', { value: 'dunnett', selected: s.multiTest === 'dunnett' }, 'One-way ANOVA + Dunnett (each vs control)')),
      h('label', null, 'Technical replicate SD: note above'), num('sdThreshold', '0.05'),
      h('label', null, 'Technical replicate SD: flag above'), num('sdMax', '0.05'),
      h('label', null, 'Late amplification: Ct above'), num('ctMax', '0.5'),
    ),
  );
}

function qcCard(res: QpcrResult) {
  const items = res.qc.map((f) => {
    let fix: Node | null = null;
    if (f.message.startsWith('Technical replicate SD')) {
      // Offer to exclude the replicate furthest from the median of the others.
      const ws = st.wells.filter((w) => !w.excluded && w.sample === f.sample && w.gene === f.gene && Number.isFinite(w.ct));
      const med = [...ws].map((w) => w.ct).sort((a, b) => a - b)[Math.floor(ws.length / 2)];
      const out = ws.reduce((a, b) => (Math.abs(b.ct - med) > Math.abs(a.ct - med) ? b : a), ws[0]);
      if (out && ws.length >= 3) fix = h('button', { class: 'small', onclick: () => update(() => (out.excluded = true)) }, `Exclude ${out.pos || '#' + (out.id + 1)} (Ct ${out.ct.toFixed(2)})`);
    }
    return h('li', { class: `flag ${f.level}` }, h('b', null, [f.sample, f.gene].filter(Boolean).join(' / ') || 'General'), f.wells.length ? ` [${f.wells.join(', ')}] ` : ' ', f.message, ' ', fix);
  });
  const excluded = st.wells.filter((w) => w.excluded);
  return card(
    '5. Quality control',
    items.length ? h('ul', { class: 'flags' }, items) : h('p', { class: 'ok' }, 'No QC flags.'),
    excluded.length ? h('p', { class: 'hint' }, `Excluded wells: ${excluded.map((w) => w.pos || '#' + (w.id + 1)).join(', ')}. `, h('button', { class: 'small', onclick: () => update(() => st.wells.forEach((w) => (w.excluded = false))) }, 'Include all again')) : null,
    res.notes.map((n) => h('p', { class: 'flag warn' }, n)),
  );
}

function perSampleCsv(res: QpcrResult) {
  return toCsv(
    ['sample', 'group', 'gene', 'ct_target_mean', 'ct_target_sd', 'ct_reference_mean', 'dCt', 'ddCt', 'fold_change_livak', 'efficiency_weighted_dCt', 'log2_ratio_pfaffl', 'ratio_pfaffl'],
    res.perSample.map((r) => [r.sample, r.group, r.gene, r.ctTarget, r.ctTargetSd, r.refCt, r.dCt, r.ddCt, r.fc, r.dCqE, r.log2Ratio, r.ratio]),
  );
}
function groupCsv(res: QpcrResult) {
  return toCsv(
    ['gene', 'group', 'n', 'dCt_mean', 'dCt_sd', 'dCt_sem', 'ddCt_mean', 'ddCt_sd', 'ddCt_sem', 'fc_mean', 'fc_sd', 'fc_sem', 'fc_of_mean_ddCt', 'pfaffl_ratio_mean', 'pfaffl_ratio_sd', 'pfaffl_ratio_sem'],
    res.groupResults.map((g) => [g.gene, g.group, g.dCt.n, g.dCt.mean, g.dCt.sd, g.dCt.sem, g.ddCt.mean, g.ddCt.sd, g.ddCt.sem, g.fc.mean, g.fc.sd, g.fc.sem, g.fcOfMean, g.ratio?.mean, g.ratio?.sd, g.ratio?.sem]),
  );
}
function statsCsv(res: QpcrResult) {
  return toCsv(['gene', 'test', 'comparison', 'difference_in_mean_dCt', 'fold_change_2^-diff', 'statistic', 'df', 'p'], res.stats.map((s) => [s.gene, s.test, s.comparison, s.diff, Number.isFinite(s.diff) ? Math.pow(2, -s.diff) : null, s.stat, String(s.df), s.p]));
}

function resultsCard(res: QpcrResult) {
  if (!res.perSample.length) return card('6. Results', h('p', { class: 'hint' }, 'Results appear once samples have groups, and at least one reference and one target gene are defined.'));
  const showPf = res.perSample.some((r) => Number.isFinite(r.ratio));
  const perRows = res.perSample.map((r) => [r.sample, r.group, h('i', null, r.gene), fmt(r.ctTarget, 2), fmt(r.refCt, 2), fmt(r.dCt), fmt(r.ddCt), fmt(r.fc), showPf ? fmt(r.ratio) : null].filter((x) => x !== null));
  const grpRows = res.groupResults.map((g) => [h('i', null, g.gene), g.group, g.dCt.n, `${fmt(g.dCt.mean)} ± ${fmt(g.dCt.sd)}`, `${fmt(g.ddCt.mean)} ± ${fmt(g.ddCt.sd)}`, `${fmt(g.fc.mean)} ± ${fmt(g.fc.sd)} (SEM ${fmt(g.fc.sem)})`, fmt(g.fcOfMean), showPf ? (g.ratio ? `${fmt(g.ratio.mean)} ± ${fmt(g.ratio.sd)}` : '') : null].filter((x) => x !== null));
  const statRows = res.stats.map((s) => [h('i', null, s.gene), s.test, s.comparison, fmt(s.diff), Number.isFinite(s.diff) ? fmt(Math.pow(2, -s.diff)) : '', fmt(s.stat), typeof s.df === 'number' ? fmt(s.df, 2) : s.df, h('b', null, formatP(s.p))]);
  return card(
    '6. Results',
    h('h3', null, 'Per biological replicate'),
    h('div', { class: 'scroll short', 'data-scroll': 'per' }, table(['Sample', 'Group', 'Gene', 'Ct target', 'Ct ref', 'ΔCt', 'ΔΔCt', '2^-ΔΔCt', ...(showPf ? ['Pfaffl ratio'] : [])], perRows)),
    h('h3', null, 'Per group (mean ± SD)'),
    table(['Gene', 'Group', 'n', 'ΔCt', 'ΔΔCt', 'Fold change (mean of replicates)', '2^-(mean ΔΔCt)', ...(showPf ? ['Pfaffl ratio'] : [])], grpRows),
    h('h3', null, 'Statistics (on ΔCt values, not on fold change)'),
    statRows.length ? table(['Gene', 'Test', 'Comparison', 'Δ mean ΔCt', 'Fold change', 'Statistic', 'df', 'p'], statRows) : h('p', { class: 'hint' }, 'No tests: at least two groups with two or more biological replicates each are needed.'),
    h('p', { class: 'hint' }, 'Δ mean ΔCt is the treatment mean minus the control mean (negative = higher expression). Welch: t statistic; Tukey: studentized range q; Dunnett: t statistic. p-values are not corrected across genes.'),
    h('div', { class: 'row' },
      h('button', { onclick: () => download('qpcr_per_sample.csv', perSampleCsv(res), 'text/csv') }, 'Per-sample CSV'),
      h('button', { onclick: () => download('qpcr_group_summary.csv', groupCsv(res), 'text/csv') }, 'Group summary CSV'),
      h('button', { onclick: () => download('qpcr_statistics.csv', statsCsv(res), 'text/csv') }, 'Statistics CSV'),
    ),
  );
}

function figureSvg(res: QpcrResult): string {
  const pf = st.settings.method === 'pfaffl';
  const series: BarSeries[] = res.targets
    .filter((g) => !st.plot.hidden.includes(g))
    .map((gene) => ({
      gene,
      groups: res.groups.map((group) => {
        const rows = res.perSample.filter((r) => r.gene === gene && r.group === group && Number.isFinite(r.ddCt));
        const usePf = pf && rows.length > 0 && rows.every((r) => Number.isFinite(r.ratio));
        const values = rows.map((r) => (st.plot.log2 ? (usePf ? r.log2Ratio : -r.ddCt) : usePf ? r.ratio : r.fc));
        const stat = res.stats.find((s) => s.gene === gene && s.groupA === group && s.groupB === res.groups[0]);
        return { group, values, p: stat?.p };
      }),
    }));
  const yLabel = st.plot.log2 ? (pf ? 'log2 ratio (Pfaffl)' : 'log2 fold change (−ΔΔCt)') : pf ? 'Relative expression (Pfaffl ratio)' : 'Fold change (2^−ΔΔCt)';
  return barPlot(series, res.groups, { log2: st.plot.log2, error: st.plot.error, showStars: st.plot.stars, showPoints: st.plot.points, title: st.plot.title, yLabel });
}

function figureCard(res: QpcrResult) {
  if (!res.perSample.length) return null;
  const svg = figureSvg(res);
  const o = st.plot;
  const chk = (k: 'log2' | 'stars' | 'points', label: string) => h('label', null, h('input', { type: 'checkbox', checked: o[k], onchange: (e: Event) => update(() => (o[k] = (e.target as HTMLInputElement).checked)) }), ' ', label);
  return card(
    '7. Figure',
    h('div', { class: 'row' },
      chk('log2', 'log2 scale'),
      chk('points', 'show replicates'),
      chk('stars', 'significance vs control'),
      h('label', null, 'Error bars ', h('select', { onchange: (e: Event) => update(() => (o.error = (e.target as HTMLSelectElement).value as 'sd' | 'sem')) }, h('option', { value: 'sd', selected: o.error === 'sd' }, 'SD'), h('option', { value: 'sem', selected: o.error === 'sem' }, 'SEM'))),
      h('input', { placeholder: 'Figure title (optional)', value: o.title, 'data-k': 'qtitle', onchange: (e: Event) => update(() => (o.title = (e.target as HTMLInputElement).value)) }),
    ),
    h('div', { class: 'row' }, h('span', null, 'Genes:'), res.targets.map((g) => h('label', null, h('input', { type: 'checkbox', checked: !o.hidden.includes(g), onchange: (e: Event) => update(() => (o.hidden = (e.target as HTMLInputElement).checked ? o.hidden.filter((x) => x !== g) : [...o.hidden, g])) }), ' ', h('i', null, g)))),
    h('div', { class: 'figure', html: svg }),
    h('p', { class: 'hint' }, 'Stars: * p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001, ns not significant; from the test chosen in step 4 (Tukey: the comparison with the control group).'),
    h('div', { class: 'row' }, h('button', { onclick: () => downloadSvg('qpcr_figure.svg', svg) }, 'Download SVG'), h('button', { onclick: () => downloadPng('qpcr_figure.png', svg) }, 'Download PNG (300 dpi)')),
  );
}
