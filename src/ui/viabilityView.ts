import absExample from '../../examples/viability_570nm_synthetic.tsv?raw';
import refExample from '../../examples/viability_650nm_reference_synthetic.tsv?raw';
import layoutExample from '../../examples/viability_layout_synthetic.json';
import { download, downloadPng, downloadSvg, toCsv } from '../lib/export';
import { parsePlateMatrix } from '../lib/parse';
import { doseResponsePlot } from '../lib/plot';
import { colorFor } from '../lib/plot';
import {
  analyzeViability,
  emptyLayout,
  fmtConc,
  layoutFromJson,
  layoutToJson,
  PLATE_COLS,
  PLATE_ROWS,
  ROW_LABELS,
  type Role,
  type ViabilityResult,
  type WellAssign,
} from '../lib/viability';
import { card, fmt, h, readFileText, rerender, table } from './dom';

type Direction = 'right' | 'left' | 'down' | 'up';

interface State {
  layout: WellAssign[][];
  unit: string;
  sel: Set<string>; // "r,c"
  absText: string;
  refText: string;
  abs: number[][] | null;
  ref: number[][] | null;
  absError: string;
  refError: string;
  fixTop: boolean;
  fixBottom: boolean;
  individual: boolean;
  title: string;
  compound: string;
  conc: string;
  topConc: string;
  factor: string;
  dir: Direction;
  example: boolean;
}

const st: State = {
  layout: emptyLayout(),
  unit: 'µM',
  sel: new Set(),
  absText: '',
  refText: '',
  abs: null,
  ref: null,
  absError: '',
  refError: '',
  fixTop: false,
  fixBottom: false,
  individual: false,
  title: '',
  compound: 'Compound 1',
  conc: '',
  topConc: '100',
  factor: '3',
  dir: 'right',
  example: false,
};

let root: HTMLElement;
let anchor: [number, number] | null = null;
let dragging = false;

export function mountViability(el: HTMLElement) {
  root = el;
  document.addEventListener('mouseup', () => (dragging = false));
  render();
}

function render() {
  rerender(root, build);
}
function update(fn: () => void) {
  fn();
  setTimeout(render, 0);
}

function readAbs() {
  st.absError = st.refError = '';
  st.abs = st.ref = null;
  if (st.absText.trim())
    try {
      st.abs = parsePlateMatrix(st.absText);
    } catch (e) {
      st.absError = (e as Error).message;
    }
  if (st.refText.trim())
    try {
      st.ref = parsePlateMatrix(st.refText);
    } catch (e) {
      st.refError = (e as Error).message;
    }
}

function loadExample() {
  const l = layoutFromJson(layoutExample);
  st.layout = l.layout;
  st.unit = l.unit;
  st.absText = absExample;
  st.refText = refExample;
  st.example = true;
  st.sel.clear();
  readAbs();
  render();
}

function build(): Node {
  const res = st.abs ? analyzeViability(st.layout, st.abs, st.ref, { top: st.fixTop ? 100 : undefined, bottom: st.fixBottom ? 0 : undefined }) : null;
  return h('div', null, introCard(), plateCard(res), absCard(), res ? [resultsCard(res), figureCard(res)] : null);
}

function introCard() {
  return card(
    'MTT / WST-1 viability and IC50',
    h('p', { class: 'hint' }, 'Define the plate layout, paste the absorbance values (and optionally a reference wavelength), and the calculator subtracts the blank, expresses each well as % of the untreated control and fits a four-parameter logistic curve per compound. One 96-well plate per analysis.'),
    h('div', { class: 'row' }, h('button', { onclick: loadExample }, 'Load synthetic example plate'), st.example ? h('span', { class: 'muted' }, 'Example loaded: simulated data for demonstration only.') : null),
  );
}

const key = (r: number, c: number) => `${r},${c}`;

function selectRect(a: [number, number], b: [number, number], add = false) {
  if (!add) st.sel.clear();
  for (let r = Math.min(a[0], b[0]); r <= Math.max(a[0], b[0]); r++) for (let c = Math.min(a[1], b[1]); c <= Math.max(a[1], b[1]); c++) st.sel.add(key(r, c));
}

function paintSelection() {
  root.querySelectorAll<HTMLElement>('.well').forEach((el) => el.classList.toggle('sel', st.sel.has(el.dataset.rc!)));
  const n = root.querySelector('#selcount');
  if (n) n.textContent = `${st.sel.size} well(s) selected`;
}

function selected(): [number, number][] {
  return [...st.sel].map((s) => s.split(',').map(Number) as [number, number]);
}

function assign(role: Role) {
  const conc = parseFloat(st.conc);
  if (role === 'compound' && (!st.compound.trim() || !(conc > 0))) return alert('Enter a compound name and a positive concentration, or use the dilution series.');
  update(() => {
    for (const [r, c] of selected()) st.layout[r][c] = { role, compound: role === 'compound' ? st.compound.trim() : '', conc: role === 'compound' ? conc : NaN, excluded: st.layout[r][c].excluded && role !== 'empty' };
  });
}

function fillSeries() {
  const top = parseFloat(st.topConc);
  const factor = parseFloat(st.factor);
  if (!st.compound.trim() || !(top > 0) || !(factor > 1)) return alert('Enter a compound name, a positive top concentration and a dilution factor above 1.');
  const cells = selected();
  if (!cells.length) return alert('Select wells first.');
  const horizontal = st.dir === 'right' || st.dir === 'left';
  const steps = [...new Set(cells.map(([r, c]) => (horizontal ? c : r)))].sort((a, b) => a - b);
  if (st.dir === 'left' || st.dir === 'up') steps.reverse();
  update(() => {
    for (const [r, c] of cells) {
      const k = steps.indexOf(horizontal ? c : r);
      st.layout[r][c] = { role: 'compound', compound: st.compound.trim(), conc: top / Math.pow(factor, k), excluded: false };
    }
  });
}

function plateCard(res: ViabilityResult | null) {
  const wellVal = (r: number, c: number) => res?.wells[r * PLATE_COLS + c];
  const compounds = [...new Set(st.layout.flat().filter((w) => w.role === 'compound').map((w) => w.compound))];
  const grid = h('div', { class: 'plate' });
  grid.append(h('div', { class: 'hdr corner', title: 'Select all', onclick: () => (selectRect([0, 0], [PLATE_ROWS - 1, PLATE_COLS - 1]), paintSelection()) }, '⌗'));
  for (let c = 0; c < PLATE_COLS; c++) grid.append(h('div', { class: 'hdr', onclick: (e: MouseEvent) => (selectRect([0, c], [PLATE_ROWS - 1, c], e.shiftKey || e.ctrlKey || e.metaKey), paintSelection()) }, String(c + 1)));
  for (let r = 0; r < PLATE_ROWS; r++) {
    grid.append(h('div', { class: 'hdr', onclick: (e: MouseEvent) => (selectRect([r, 0], [r, PLATE_COLS - 1], e.shiftKey || e.ctrlKey || e.metaKey), paintSelection()) }, ROW_LABELS[r]));
    for (let c = 0; c < PLATE_COLS; c++) {
      const a = st.layout[r][c];
      const v = wellVal(r, c);
      const ci = compounds.indexOf(a.compound);
      const style = a.role === 'compound' ? `--wc:${colorFor(ci)}` : '';
      const label = a.role === 'blank' ? 'Blank' : a.role === 'control' ? 'Ctrl' : a.role === 'compound' ? fmtConc(a.conc) : '';
      grid.append(
        h('div', {
          class: `well ${a.role}${a.excluded ? ' excl' : ''}${st.sel.has(key(r, c)) ? ' sel' : ''}`,
          style,
          'data-rc': key(r, c),
          title: `${ROW_LABELS[r]}${c + 1}: ${a.role}${a.role === 'compound' ? ` ${a.compound} ${a.conc} ${st.unit}` : ''}${a.excluded ? ' (excluded)' : ''}${v && Number.isFinite(v.raw) ? `\nA = ${v.raw}${st.ref ? `, corrected ${fmt(v.corrected, 4)}` : ''}` : ''}${v && Number.isFinite(v.viability) ? `\n${v.viability.toFixed(1)} % viability` : ''}`,
          onmousedown: (e: MouseEvent) => {
            e.preventDefault();
            dragging = true;
            if (e.shiftKey && anchor) selectRect(anchor, [r, c]);
            else if (e.ctrlKey || e.metaKey) {
              anchor = [r, c];
              const k = key(r, c);
              if (st.sel.has(k)) st.sel.delete(k);
              else st.sel.add(k);
            } else {
              anchor = [r, c];
              selectRect([r, c], [r, c]);
            }
            paintSelection();
          },
          onmouseenter: () => {
            if (dragging && anchor) {
              selectRect(anchor, [r, c]);
              paintSelection();
            }
          },
        },
          h('span', { class: 'lbl' }, label),
          v && Number.isFinite(v.raw) ? h('span', { class: 'od' }, v.raw.toFixed(3)) : null,
        ),
      );
    }
  }
  return card(
    '1. Plate layout',
    h('p', { class: 'hint' }, 'Click or drag to select wells (Shift extends, Ctrl/Cmd adds; click a row or column header to select it), then assign a role. Excluded wells keep their role but are left out of every calculation.'),
    h('div', { class: 'plate-wrap' },
      grid,
      h('div', { class: 'assign' },
        h('div', { id: 'selcount', class: 'muted' }, `${st.sel.size} well(s) selected`),
        compounds.length ? h('div', { class: 'row' }, compounds.map((c, i) => h('span', null, h('span', { class: 'swatch', style: `background:${colorFor(i)}` }), ' ', c))) : null,
        h('div', { class: 'row' }, h('button', { onclick: () => assign('blank') }, 'Blank (medium only)'), h('button', { onclick: () => assign('control') }, 'Untreated control')),
        h('div', { class: 'row' }, h('button', { onclick: () => assign('empty') }, 'Clear'), h('button', { onclick: () => update(() => selected().forEach(([r, c]) => (st.layout[r][c].excluded = !st.layout[r][c].excluded))) }, 'Toggle exclude')),
        h('fieldset', null, h('legend', null, 'Compound'),
          h('label', null, 'Name ', h('input', { value: st.compound, list: 'cmpdlist', onchange: (e: Event) => (st.compound = (e.target as HTMLInputElement).value) })),
          h('datalist', { id: 'cmpdlist' }, compounds.map((c) => h('option', { value: c }))),
          h('label', null, 'Unit ', h('select', { onchange: (e: Event) => update(() => (st.unit = (e.target as HTMLSelectElement).value)) }, ['nM', 'µM', 'mM', 'ng/mL', 'µg/mL', 'mg/mL', '%'].map((u) => h('option', { value: u, selected: u === st.unit }, u)))),
          h('div', { class: 'row' }, h('label', null, 'Concentration ', h('input', { class: 'short', type: 'number', value: st.conc, onchange: (e: Event) => (st.conc = (e.target as HTMLInputElement).value) })), h('button', { onclick: () => assign('compound') }, 'Set')),
          h('div', { class: 'row' },
            h('label', null, 'Top ', h('input', { class: 'short', type: 'number', value: st.topConc, onchange: (e: Event) => (st.topConc = (e.target as HTMLInputElement).value) })),
            h('label', null, '÷ ', h('input', { class: 'short', type: 'number', value: st.factor, onchange: (e: Event) => (st.factor = (e.target as HTMLInputElement).value) })),
          ),
          h('div', { class: 'row' },
            h('select', { onchange: (e: Event) => (st.dir = (e.target as HTMLSelectElement).value as Direction) }, ([['right', 'decreasing left → right'], ['left', 'decreasing right → left'], ['down', 'decreasing top → bottom'], ['up', 'decreasing bottom → top']] as const).map(([v, t]) => h('option', { value: v, selected: st.dir === v }, t))),
            h('button', { onclick: fillSeries }, 'Fill dilution series'),
          ),
        ),
        h('div', { class: 'row' },
          h('button', { onclick: () => download('plate_layout.json', JSON.stringify(layoutToJson(st.layout, st.unit), null, 1), 'application/json') }, 'Save layout'),
          h('label', { class: 'btn' }, 'Load layout', h('input', {
            type: 'file',
            accept: '.json',
            onchange: async (e: Event) => {
              const f = (e.target as HTMLInputElement).files?.[0];
              if (!f) return;
              try {
                const l = layoutFromJson(JSON.parse(await readFileText(f)));
                update(() => ((st.layout = l.layout), (st.unit = l.unit)));
              } catch (err) {
                alert((err as Error).message);
              }
            },
          })),
          h('button', { onclick: () => update(() => ((st.layout = emptyLayout()), st.sel.clear())) }, 'Clear plate'),
        ),
      ),
    ),
  );
}

function absCard() {
  return card(
    '2. Absorbance',
    h('p', { class: 'hint' }, 'Paste the 8 × 12 block of values from the plate reader or a spreadsheet (row letters and a 1 to 12 header row are optional; decimal commas are accepted). MTT is typically read at 570 nm, WST-1 at 440 to 450 nm.'),
    h('div', { class: 'cols' },
      h('div', null, h('h3', null, 'Measurement wavelength'), h('textarea', { rows: 9, value: st.absText, onchange: (e: Event) => update(() => ((st.absText = (e.target as HTMLTextAreaElement).value), readAbs())) }), st.absError ? h('p', { class: 'flag error' }, st.absError) : st.abs ? h('p', { class: 'ok' }, 'Read 96 values.') : null),
      h('div', null, h('h3', null, 'Reference wavelength (optional)'), h('textarea', { rows: 9, value: st.refText, placeholder: 'e.g. 630 or 650 nm; subtracted well by well', onchange: (e: Event) => update(() => ((st.refText = (e.target as HTMLTextAreaElement).value), readAbs())) }), st.refError ? h('p', { class: 'flag error' }, st.refError) : st.ref ? h('p', { class: 'ok' }, 'Reference values will be subtracted.') : null),
    ),
  );
}

function wellsCsv(res: ViabilityResult) {
  return toCsv(['well', 'role', 'compound', `conc_${st.unit}`, 'excluded', 'absorbance', 'reference_corrected', 'viability_percent'], res.wells.filter((w) => w.role !== 'empty').map((w) => [w.pos, w.role, w.compound, Number.isFinite(w.conc) ? w.conc : null, w.excluded ? 'yes' : 'no', w.raw, w.corrected, w.viability]));
}
function replicateCsv(res: ViabilityResult) {
  return toCsv(['compound', `conc_${st.unit}`, 'n', 'viability_mean', 'viability_sd'], res.compounds.flatMap((c) => c.points.map((p) => [c.name, p.conc, p.values.length, p.mean, p.sd])));
}
function fitCsv(res: ViabilityResult) {
  return toCsv(
    ['compound', 'n', 'top', 'top_se', 'bottom', 'bottom_se', 'hill', 'hill_se', 'logIC50', 'logIC50_se', `IC50_${st.unit}`, 'IC50_CI95_low', 'IC50_CI95_high', `absolute_IC50_${st.unit}`, 'R2', 'top_fixed', 'bottom_fixed', 'flags'],
    res.compounds.map((c) => {
      const f = c.fit;
      return [c.name, f?.n, f?.params.top, f?.se.top, f?.params.bottom, f?.se.bottom, f?.params.hill, f?.se.hill, f?.params.logIC50, f?.se.logIC50, f?.ic50, f?.ic50CI[0], f?.ic50CI[1], f?.absIC50, f?.r2, st.fixTop ? 'yes' : 'no', st.fixBottom ? 'yes' : 'no', c.flags.join(' | ')];
    }),
  );
}

function resultsCard(res: ViabilityResult) {
  const fitRows = res.compounds.map((c, i) => {
    const f = c.fit;
    return [
      h('span', { class: 'swatch', style: `background:${colorFor(i)}` }), c.name,
      f ? f.n : '',
      f ? fmt(f.params.top, 1) : '',
      f ? fmt(f.params.bottom, 1) : '',
      f ? fmt(f.params.hill, 2) : '',
      f ? h('b', null, fmtConc(f.ic50)) : '',
      f ? (f.ic50CI[0] > 0 && Number.isFinite(f.ic50CI[1]) ? `${fmtConc(f.ic50CI[0])} to ${fmtConc(f.ic50CI[1])}` : 'not determinable') : '',
      f ? fmtConc(f.absIC50) : '',
      f ? fmt(f.r2, 4) : '',
      c.flags.length ? h('ul', { class: 'flags inline' }, c.flags.map((x) => h('li', { class: 'flag warn' }, x))) : h('span', { class: 'ok' }, 'none'),
    ];
  });
  const repRows = res.compounds.flatMap((c) => c.points.map((p) => [c.name, fmtConc(p.conc), p.values.length, fmt(p.mean, 1), fmt(p.sd, 1)]));
  return card(
    '3. Results',
    h('div', { class: 'row' },
      h('label', null, h('input', { type: 'checkbox', checked: st.fixTop, onchange: (e: Event) => update(() => (st.fixTop = (e.target as HTMLInputElement).checked)) }), ' Fix top plateau at 100 %'),
      h('label', null, h('input', { type: 'checkbox', checked: st.fixBottom, onchange: (e: Event) => update(() => (st.fixBottom = (e.target as HTMLInputElement).checked)) }), ' Fix bottom plateau at 0 %'),
    ),
    h('p', { class: 'hint' }, 'Default is the free four-parameter fit. Fix a plateau only when the data do not define it (for example, no complete kill at the highest dose) and you have a reason to assume its value.'),
    res.warnings.map((w) => h('p', { class: 'flag warn' }, w)),
    h('p', null, `Blank: mean ${fmt(res.blankMean, 4)} (SD ${fmt(res.blankSd, 4)}, n = ${res.nBlank}). Untreated control after blank subtraction: mean ${fmt(res.controlMean, 4)} (SD ${fmt(res.controlSd, 4)}, n = ${res.nControl}) = 100 %.`),
    h('h3', null, 'Dose-response fits (four-parameter logistic)'),
    table(['', 'Compound', 'n', 'Top %', 'Bottom %', 'Hill', `IC50 (${st.unit})`, '95% CI', `Abs. IC50 (${st.unit})`, 'R²', 'Flags'], fitRows),
    h('p', { class: 'hint' }, 'IC50 is the relative IC50 (midpoint between the fitted plateaus); its 95% CI is computed on the log scale and back-transformed. Absolute IC50 is the concentration where the curve crosses 50 % viability.'),
    h('details', null, h('summary', null, 'Replicate means'), table(['Compound', `Conc (${st.unit})`, 'n', 'Mean %', 'SD %'], repRows)),
    h('div', { class: 'row' },
      h('button', { onclick: () => download('viability_wells.csv', wellsCsv(res), 'text/csv') }, 'Per-well CSV'),
      h('button', { onclick: () => download('viability_replicates.csv', replicateCsv(res), 'text/csv') }, 'Replicate summary CSV'),
      h('button', { onclick: () => download('viability_fits.csv', fitCsv(res), 'text/csv') }, 'Fit parameters CSV'),
    ),
  );
}

function figureCard(res: ViabilityResult) {
  if (!res.compounds.length) return null;
  const svg = doseResponsePlot(res.compounds, { unit: st.unit, individual: st.individual, title: st.title });
  return card(
    '4. Figure',
    h('div', { class: 'row' },
      h('label', null, h('input', { type: 'checkbox', checked: st.individual, onchange: (e: Event) => update(() => (st.individual = (e.target as HTMLInputElement).checked)) }), ' show individual wells instead of mean ± SD'),
      h('input', { placeholder: 'Figure title (optional)', value: st.title, 'data-k': 'vtitle', onchange: (e: Event) => update(() => (st.title = (e.target as HTMLInputElement).value)) }),
    ),
    h('div', { class: 'figure', html: svg }),
    h('div', { class: 'row' }, h('button', { onclick: () => downloadSvg('dose_response.svg', svg) }, 'Download SVG'), h('button', { onclick: () => downloadPng('dose_response.png', svg) }, 'Download PNG (300 dpi)')),
  );
}
