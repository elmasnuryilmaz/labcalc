// MTT / WST-1 viability: reference-wavelength and blank subtraction, % viability
// relative to untreated control, replicate summaries and 4PL dose-response fits.
import { fit4PL, type FitResult, type FixedParams } from './fit4pl';
import { mean, sd } from './stats';

export type Role = 'empty' | 'blank' | 'control' | 'compound';

export interface WellAssign {
  role: Role;
  compound: string;
  conc: number;
  excluded: boolean;
}

export const PLATE_ROWS = 8;
export const PLATE_COLS = 12;
export const ROW_LABELS = 'ABCDEFGH';

export function emptyLayout(): WellAssign[][] {
  return Array.from({ length: PLATE_ROWS }, () =>
    Array.from({ length: PLATE_COLS }, () => ({ role: 'empty' as Role, compound: '', conc: NaN, excluded: false })),
  );
}

export interface DosePoint {
  conc: number;
  values: number[];
  mean: number;
  sd: number;
}

export interface CompoundResult {
  name: string;
  points: DosePoint[];
  fit: FitResult | null;
  flags: string[];
}

export interface WellValue {
  pos: string;
  row: number;
  col: number;
  role: Role;
  compound: string;
  conc: number;
  excluded: boolean;
  raw: number;
  corrected: number; // after reference-wavelength subtraction
  viability: number; // percent of untreated control
}

export interface ViabilityResult {
  blankMean: number;
  blankSd: number;
  nBlank: number;
  controlMean: number; // blank-subtracted
  controlSd: number;
  nControl: number;
  wells: WellValue[];
  compounds: CompoundResult[];
  warnings: string[];
}

export function analyzeViability(layout: WellAssign[][], abs: number[][], ref: number[][] | null, fixed: FixedParams = {}): ViabilityResult {
  const warnings: string[] = [];
  const wells: WellValue[] = [];
  for (let r = 0; r < PLATE_ROWS; r++)
    for (let c = 0; c < PLATE_COLS; c++) {
      const a = layout[r][c];
      const raw = abs[r]?.[c] ?? NaN;
      const corrected = ref ? raw - (ref[r]?.[c] ?? NaN) : raw;
      wells.push({ pos: `${ROW_LABELS[r]}${c + 1}`, row: r, col: c, role: a.role, compound: a.compound, conc: a.conc, excluded: a.excluded, raw, corrected, viability: NaN });
    }
  const use = (w: WellValue) => !w.excluded && w.role !== 'empty';
  const missing = wells.filter((w) => use(w) && !Number.isFinite(w.corrected));
  if (missing.length) warnings.push(`${missing.length} assigned well(s) have no absorbance value: ${missing.map((w) => w.pos).join(', ')}.`);

  const blanks = wells.filter((w) => use(w) && w.role === 'blank' && Number.isFinite(w.corrected)).map((w) => w.corrected);
  let blankMean = mean(blanks);
  if (blanks.length === 0) {
    warnings.push('No blank (medium only) wells: absorbance is used without blank subtraction.');
    blankMean = 0;
  }
  const controls = wells.filter((w) => use(w) && w.role === 'control' && Number.isFinite(w.corrected)).map((w) => w.corrected - blankMean);
  const controlMean = mean(controls);
  const res: ViabilityResult = { blankMean, blankSd: sd(blanks), nBlank: blanks.length, controlMean, controlSd: sd(controls), nControl: controls.length, wells, compounds: [], warnings };
  if (controls.length === 0) {
    warnings.push('No untreated control wells: % viability cannot be calculated.');
    return res;
  }
  if (controlMean <= 0) warnings.push('Mean blank-subtracted control absorbance is not positive; check blank and control wells.');
  if (controls.length >= 2 && sd(controls) / controlMean > 0.15) warnings.push(`Untreated control CV is ${((100 * sd(controls)) / controlMean).toFixed(1)}% (above 15%).`);

  for (const w of wells) if (Number.isFinite(w.corrected) && w.role !== 'empty') w.viability = (100 * (w.corrected - blankMean)) / controlMean;

  const names = [...new Set(wells.filter((w) => w.role === 'compound' && w.compound).map((w) => w.compound))];
  for (const name of names) {
    const cw = wells.filter((w) => use(w) && w.role === 'compound' && w.compound === name && Number.isFinite(w.viability));
    const concs = [...new Set(cw.map((w) => w.conc).filter((c) => Number.isFinite(c) && c > 0))].sort((a, b) => a - b);
    const points = concs.map((c) => {
      const v = cw.filter((w) => w.conc === c).map((w) => w.viability);
      return { conc: c, values: v, mean: mean(v), sd: sd(v) };
    });
    const flags: string[] = [];
    if (cw.some((w) => !(w.conc > 0))) flags.push('Wells without a positive concentration are not used in the fit.');
    const xs = points.flatMap((p) => p.values.map(() => p.conc));
    const ys = points.flatMap((p) => p.values);
    const fit = points.length >= 3 ? fit4PL(xs, ys, fixed) : null;
    if (!fit) flags.push('Not enough data to fit (needs at least 3 concentrations and one value per fitted parameter).');
    else {
      const lo = concs[0];
      const hi = concs[concs.length - 1];
      if (!fit.converged) flags.push('Fit did not converge; treat the parameters with caution.');
      if (fit.ic50 < lo || fit.ic50 > hi) flags.push(`IC50 lies outside the tested range (${fmtConc(lo)} to ${fmtConc(hi)}); it is an extrapolation.`);
      if (fit.df <= 0) flags.push('No residual degrees of freedom: confidence intervals are not available.');
      else if (fit.logIC50CI[1] - fit.logIC50CI[0] > 2) flags.push('95% CI of the IC50 spans more than two orders of magnitude.');
      if (fit.params.bottom > fit.params.top) flags.push('Response increases with concentration (bottom plateau above top plateau).');
      if (fit.params.hill > 5) flags.push('Very steep Hill slope (> 5): a plateau is probably not defined by the data. Consider fixing the bottom or top plateau.');
      if (!Number.isFinite(fit.absIC50)) flags.push('Fitted curve does not cross 50% viability; the IC50 reported is the relative (midpoint) IC50.');
    }
    res.compounds.push({ name, points, fit, flags });
  }
  return res;
}

export function fmtConc(c: number): string {
  if (!Number.isFinite(c)) return 'n/a';
  if (c >= 1000 || c < 0.01) return c.toExponential(2);
  return Number(c.toPrecision(3)).toString();
}

export interface LayoutFile {
  format: 'labcalc-plate-layout';
  version: 1;
  unit: string;
  note?: string;
  wells: { pos: string; role: Role; compound?: string; conc?: number; excluded?: boolean }[];
}

/** Read a plate layout saved with layoutToJson (or the synthetic example layout). */
export function layoutFromJson(obj: unknown): { layout: WellAssign[][]; unit: string } {
  const f = obj as LayoutFile;
  if (!f || f.format !== 'labcalc-plate-layout' || !Array.isArray(f.wells)) throw new Error('Not a labcalc plate layout file.');
  const layout = emptyLayout();
  for (const w of f.wells) {
    const m = String(w.pos).toUpperCase().match(/^([A-H])(\d{1,2})$/);
    if (!m) continue;
    const r = ROW_LABELS.indexOf(m[1]);
    const c = Number(m[2]) - 1;
    if (c < 0 || c >= PLATE_COLS) continue;
    const role: Role = ['blank', 'control', 'compound'].includes(w.role) ? w.role : 'empty';
    layout[r][c] = { role, compound: role === 'compound' ? String(w.compound ?? '') : '', conc: role === 'compound' ? Number(w.conc) : NaN, excluded: !!w.excluded };
  }
  return { layout, unit: f.unit || 'µM' };
}

export function layoutToJson(layout: WellAssign[][], unit: string): LayoutFile {
  const wells: LayoutFile['wells'] = [];
  layout.forEach((row, r) =>
    row.forEach((a, c) => {
      if (a.role === 'empty' && !a.excluded) return;
      const w: LayoutFile['wells'][number] = { pos: `${ROW_LABELS[r]}${c + 1}`, role: a.role };
      if (a.role === 'compound') {
        w.compound = a.compound;
        w.conc = a.conc;
      }
      if (a.excluded) w.excluded = true;
      wells.push(w);
    }),
  );
  return { format: 'labcalc-plate-layout', version: 1, unit, wells };
}
