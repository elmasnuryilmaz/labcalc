import { defineConfig, type Plugin } from 'vite';

// Strict Content-Security-Policy for the production build: no network requests
// of any kind (connect-src 'none'), no external scripts, styles, fonts or images.
const csp = [
  "default-src 'none'",
  "script-src 'self'",
  "style-src 'self' 'unsafe-inline'",
  "img-src 'self' data: blob:",
  "connect-src 'none'",
  "font-src 'none'",
  "object-src 'none'",
  "base-uri 'none'",
  "form-action 'none'",
].join('; ');

const cspPlugin = (): Plugin => ({
  name: 'labcalc-csp',
  apply: 'build',
  transformIndexHtml: (html) => html.replace('<meta charset="UTF-8" />', `<meta charset="UTF-8" />\n    <meta http-equiv="Content-Security-Policy" content="${csp}" />`),
});

export default defineConfig({
  base: './',
  plugins: [cspPlugin()],
  build: { target: 'es2020', assetsInlineLimit: 0 },
  test: { include: ['tests/**/*.test.ts'] },
});
