import { describe, expect, it } from 'vitest';
import { fmt } from '../src/ui/dom';

describe('fmt', () => {
  it('shows floating-point residue as zero', () => {
    expect(fmt(-8.881784197001252e-16)).toBe('0.000');
    expect(fmt(4.4e-13, 2)).toBe('0.00');
  });

  it('keeps small and large real values in exponential form', () => {
    expect(fmt(0.00042)).toBe('4.20e-4');
    expect(fmt(123456)).toBe('1.23e+5');
  });

  it('formats ordinary values and drops non-finite ones', () => {
    expect(fmt(1.9234)).toBe('1.923');
    expect(fmt(-0.565, 2)).toBe('-0.56');
    expect(fmt(NaN)).toBe('');
  });
});
