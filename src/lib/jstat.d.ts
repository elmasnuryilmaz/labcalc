// Minimal type declarations for the parts of jStat used by labcalc.
declare module 'jstat' {
  interface JStatStatic {
    ibeta(x: number, a: number, b: number): number;
    gammaln(x: number): number;
    normal: { cdf(x: number, mean: number, sd: number): number };
    studentt: { cdf(x: number, df: number): number; inv(p: number, df: number): number };
    tukey: { cdf(q: number, nmeans: number, df: number): number };
  }
  const jStat: JStatStatic;
  export default jStat;
}
