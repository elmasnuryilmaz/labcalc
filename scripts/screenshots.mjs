// Drive the built app in headless Chromium: load the synthetic examples, check that
// the page makes no requests outside its own origin and has no console errors,
// exercise the PNG export, and save the README screenshots.
// Usage: npm run build && npx vite preview --port 4173 & node scripts/screenshots.mjs
import { createRequire } from 'node:module';

// Uses a globally or locally installed Playwright (set NODE_PATH to the global node_modules if needed).
const { chromium } = createRequire(import.meta.url)('playwright');

const base = process.env.LABCALC_URL ?? 'http://localhost:4173/';
const out = 'docs/screenshots';
const browser = await chromium.launch();
const page = await browser.newPage({ viewport: { width: 1280, height: 900 }, deviceScaleFactor: 1, acceptDownloads: true });
const errors = [];
const foreign = [];
page.on('console', (m) => m.type() === 'error' && errors.push(m.text()));
page.on('pageerror', (e) => errors.push(String(e)));
page.on('request', (r) => !r.url().startsWith(base) && !r.url().startsWith('blob:') && !r.url().startsWith('data:') && foreign.push(r.url()));

const panel = () => page.locator('[role=tabpanel]:not([hidden])');
await page.goto(base + '#qpcr');
await page.getByRole('button', { name: /synthetic LC480 example/ }).click();
await panel().locator('.figure svg').waitFor();
const card = (title) => panel().locator('section.card', { has: page.locator('h2', { hasText: title }) });
await card('5. Quality control').screenshot({ path: `${out}/qpcr-qc.png` });
await card('6. Results').screenshot({ path: `${out}/qpcr-results.png` });
await card('7. Figure').locator('.figure').screenshot({ path: `${out}/qpcr-figure.png` });
const dl = page.waitForEvent('download');
await card('7. Figure').getByRole('button', { name: /PNG/ }).click();
const png = await dl;
console.log('PNG export:', png.suggestedFilename(), (await png.path()) ? 'ok' : 'failed');

await page.goto(base + '#viability');
await page.getByRole('button', { name: /synthetic example plate/ }).click();
await panel().locator('.figure svg').waitFor();
await card('1. Plate layout').screenshot({ path: `${out}/viability-plate.png` });
await card('3. Results').screenshot({ path: `${out}/viability-results.png` });
await card('4. Figure').locator('.figure').screenshot({ path: `${out}/viability-figure.png` });
const text = await card('3. Results').innerText();
console.log(text.split('\n').filter((l) => /Cmpd|Blank/.test(l)).slice(0, 8).join('\n'));

console.log('console errors:', errors.length ? errors : 'none');
console.log('requests outside the page origin:', foreign.length ? foreign : 'none');
await browser.close();
if (errors.length || foreign.length) process.exit(1);
