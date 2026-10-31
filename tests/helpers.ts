import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { expect } from 'vitest';

export const root = resolve(__dirname, '..');
export const readText = (p: string) => readFileSync(resolve(root, p), 'utf-8');
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export const ref: any = JSON.parse(readText('tests/fixtures/reference.json'));

/** |actual - expected| <= atol + rtol * |expected|; null in the fixture means NaN. */
export function close(actual: number, expected: number | null, rtol = 1e-6, atol = 1e-9, label = '') {
  if (expected === null) {
    expect(Number.isFinite(actual), `${label}: expected NaN/inf, got ${actual}`).toBe(false);
    return;
  }
  const tol = atol + rtol * Math.abs(expected);
  expect(Math.abs(actual - expected), `${label}: actual ${actual}, expected ${expected}`).toBeLessThanOrEqual(tol);
}
