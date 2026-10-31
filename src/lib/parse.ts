// Parsers for qPCR exports (generic CSV/TSV and Roche LightCycler 480 text export)
// and for 96-well plate absorbance matrices pasted from a plate reader or spreadsheet.

export interface Well {
  id: number; // running index, unique within a dataset
  pos: string; // well position, e.g. "A1" (may be empty for generic tables)
  name: string; // sample name as it appears in the file
  ct: number; // NaN when undetermined / missing
  ctRaw: string; // original text of the Ct cell
  sample: string;
  gene: string;
  group: string; // only used when the file provides a group column
  excluded: boolean;
}

/** Decode file bytes, handling UTF-8 and UTF-16 (LE/BE) byte-order marks. */
export function decodeBytes(buf: ArrayBuffer): string {
  const b = new Uint8Array(buf);
  if (b.length >= 2 && b[0] === 0xff && b[1] === 0xfe) return new TextDecoder('utf-16le').decode(b.subarray(2));
  if (b.length >= 2 && b[0] === 0xfe && b[1] === 0xff) return new TextDecoder('utf-16be').decode(b.subarray(2));
  // UTF-16 without BOM: many zero bytes in odd positions.
  if (b.length > 4 && b[1] === 0 && b[3] === 0) return new TextDecoder('utf-16le').decode(b);
  return new TextDecoder('utf-8').decode(b).replace(/^﻿/, '');
}

/** Parse a number, accepting a decimal comma. Returns NaN for anything non-numeric. */
export function parseNumber(s: string | undefined): number {
  if (s === undefined) return NaN;
  let t = s.trim().replace(/^"|"$/g, '').trim();
  if (t === '') return NaN;
  if (/^-?\d+,\d+$/.test(t)) t = t.replace(',', '.');
  if (!/^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(t)) return NaN;
  return Number(t);
}

function splitLine(line: string, delim: string): string[] {
  if (delim === '\t') return line.split('\t').map((c) => c.trim());
  // Minimal quoted-field CSV splitting.
  const out: string[] = [];
  let cur = '';
  let q = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (q) {
      if (ch === '"' && line[i + 1] === '"') {
        cur += '"';
        i++;
      } else if (ch === '"') q = false;
      else cur += ch;
    } else if (ch === '"') q = true;
    else if (ch === delim) {
      out.push(cur.trim());
      cur = '';
    } else cur += ch;
  }
  out.push(cur.trim());
  return out;
}

export function detectDelimiter(line: string): string {
  const counts: [string, number][] = ['\t', ';', ','].map((d) => [d, line.split(d).length - 1]);
  counts.sort((a, b) => b[1] - a[1]);
  return counts[0][1] > 0 ? counts[0][0] : '\t';
}

export function splitRows(text: string): string[] {
  return text.replace(/\r\n?/g, '\n').split('\n');
}

const SYN = {
  pos: ['well', 'pos', 'position', 'well position'],
  name: ['name', 'sample name', 'sample'],
  sample: ['sample', 'sample name', 'biological sample', 'bio sample'],
  gene: ['gene', 'target', 'target name', 'detector', 'assay', 'primer'],
  group: ['group', 'condition', 'treatment'],
  ct: ['ct', 'cp', 'cq', 'c(t)', 'crt'],
};

function findCol(header: string[], keys: string[]): number {
  const h = header.map((x) => x.toLowerCase().replace(/^"|"$/g, '').trim());
  for (const k of keys) {
    const i = h.indexOf(k);
    if (i >= 0) return i;
  }
  return -1;
}

export interface QpcrParseResult {
  format: 'lc480' | 'table';
  wells: Well[];
  warnings: string[];
}

/**
 * Parse a qPCR export. Recognises the LightCycler 480 "Abs Quant / 2nd Derivative Max"
 * text export (optional "Experiment:" line, then a tab-separated header with
 * Include, Color, Pos, Name, Cp, Concentration, Standard, Status) and generic
 * CSV/TSV tables with a header row containing a Ct/Cp/Cq column.
 */
export function parseQpcr(text: string): QpcrParseResult {
  const lines = splitRows(text);
  const warnings: string[] = [];
  let headerIdx = -1;
  let delim = '\t';
  for (let i = 0; i < Math.min(lines.length, 50); i++) {
    const d = detectDelimiter(lines[i]);
    const cells = splitLine(lines[i], d);
    if (findCol(cells, SYN.ct) >= 0) {
      headerIdx = i;
      delim = d;
      break;
    }
  }
  if (headerIdx < 0) throw new Error('No header row with a Ct, Cp or Cq column was found.');
  const header = splitLine(lines[headerIdx], delim);
  const lower = header.map((h) => h.toLowerCase());
  const isLC = lower.includes('pos') && lower.includes('cp') && (lower.includes('include') || lower.includes('status'));

  const cPos = findCol(header, SYN.pos);
  const cCt = findCol(header, SYN.ct);
  const cInclude = findCol(header, ['include']);
  const cSample = isLC ? -1 : findCol(header, SYN.sample);
  const cName = findCol(header, SYN.name);
  const cGene = findCol(header, SYN.gene);
  const cGroup = findCol(header, SYN.group);

  const wells: Well[] = [];
  for (let i = headerIdx + 1; i < lines.length; i++) {
    const raw = lines[i];
    if (raw.trim() === '') continue;
    const c = splitLine(raw, delim);
    const get = (idx: number) => (idx >= 0 && idx < c.length ? c[idx] : '');
    const ctRaw = get(cCt);
    const name = get(cName >= 0 ? cName : cSample);
    const pos = get(cPos).toUpperCase();
    if (!pos && !name && !ctRaw) continue;
    const include = cInclude >= 0 ? !/^false$/i.test(get(cInclude)) : true;
    wells.push({
      id: wells.length,
      pos,
      name,
      ct: parseNumber(ctRaw),
      ctRaw,
      sample: cSample >= 0 ? get(cSample) : '',
      gene: cGene >= 0 ? get(cGene) : '',
      group: cGroup >= 0 ? get(cGroup) : '',
      excluded: !include,
    });
  }
  if (wells.length === 0) warnings.push('The file contains a header but no data rows.');
  const notIncluded = wells.filter((w) => w.excluded).length;
  if (notIncluded > 0) warnings.push(`${notIncluded} well(s) had Include = False in the export and start excluded.`);
  return { format: isLC ? 'lc480' : 'table', wells, warnings };
}

/** Split well names such as "Ctrl1_GAPDH" into sample and gene parts. */
export function splitName(name: string, delim: string, genePart: 'first' | 'last'): { sample: string; gene: string } {
  const d = delim === '' ? /\s+/ : delim;
  const parts = name.trim().split(d).filter((p) => p !== '');
  if (parts.length < 2) return { sample: name.trim(), gene: '' };
  if (genePart === 'last') return { gene: parts[parts.length - 1], sample: parts.slice(0, -1).join(typeof d === 'string' ? d : ' ') };
  return { gene: parts[0], sample: parts.slice(1).join(typeof d === 'string' ? d : ' ') };
}

/** Guess a group label from a sample name by removing a trailing replicate number, e.g. "Ctrl-2" -> "Ctrl". */
export function guessGroup(sample: string): string {
  const g = sample.replace(/[\s_\-.]*(?:rep[\s_\-.]*|#)?\d+$/i, '').trim();
  return g || sample;
}

export const ROWS = 'ABCDEFGHIJKLMNOP';

/** Parse "A1-B6, C3" into a list of well positions (rectangular blocks). */
export function parseWellRanges(spec: string): string[] {
  const out: string[] = [];
  for (const part of spec.split(/[,;\s]+/).filter(Boolean)) {
    const m = part.toUpperCase().match(/^([A-P])(\d{1,2})(?:[-:]([A-P])(\d{1,2}))?$/);
    if (!m) throw new Error(`Cannot read well range "${part}". Use e.g. A1-A12 or A1:B6.`);
    const r1 = ROWS.indexOf(m[1]);
    const c1 = Number(m[2]);
    const r2 = m[3] ? ROWS.indexOf(m[3]) : r1;
    const c2 = m[4] ? Number(m[4]) : c1;
    for (let r = Math.min(r1, r2); r <= Math.max(r1, r2); r++)
      for (let c = Math.min(c1, c2); c <= Math.max(c1, c2); c++) out.push(`${ROWS[r]}${c}`);
  }
  return out;
}

/** Normalise "A01" to "A1". */
export function normPos(pos: string): string {
  const m = pos.toUpperCase().match(/^([A-P])0*(\d{1,2})$/);
  return m ? `${m[1]}${Number(m[2])}` : pos.toUpperCase();
}

/**
 * Parse an 8 x 12 absorbance matrix pasted from a spreadsheet or plate reader.
 * Accepts an optional header row (1..12) and an optional row-label column (A..H).
 */
export function parsePlateMatrix(text: string, rows = 8, cols = 12): number[][] {
  const lines = splitRows(text).filter((l) => l.trim() !== '');
  const matrix: number[][] = [];
  for (const line of lines) {
    const d = detectDelimiter(line);
    let cells = d === '\t' && !line.includes('\t') ? line.trim().split(/\s+/) : splitLine(line, d);
    if (cells.length > 0 && /^[A-P]$/i.test(cells[0].trim())) cells = cells.slice(1);
    const nums = cells.map(parseNumber);
    const numeric = nums.filter((v) => Number.isFinite(v));
    // Skip a header row of column numbers 1..12.
    if (numeric.length >= cols && numeric.slice(0, cols).every((v, i) => v === i + 1) && matrix.length === 0) continue;
    if (numeric.length === 0) continue;
    matrix.push(nums.slice(0, cols).concat(new Array(Math.max(0, cols - nums.length)).fill(NaN)));
    if (matrix.length === rows) break;
  }
  if (matrix.length !== rows) throw new Error(`Expected ${rows} rows of absorbance values, found ${matrix.length}.`);
  return matrix;
}
