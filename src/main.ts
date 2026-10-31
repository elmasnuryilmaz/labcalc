import './style.css';
import { h } from './ui/dom';
import { mountAbout, mountMethods } from './ui/methodsView';
import { mountQpcr } from './ui/qpcrView';
import { mountViability } from './ui/viabilityView';

const tabs = [
  { id: 'qpcr', label: 'qPCR ΔΔCt', mount: mountQpcr },
  { id: 'viability', label: 'Cell viability (MTT / WST-1)', mount: mountViability },
  { id: 'methods', label: 'Methods', mount: mountMethods },
  { id: 'about', label: 'Privacy & about', mount: mountAbout },
];

const panels = new Map<string, HTMLElement>();
const nav = h('nav', { class: 'tabs', role: 'tablist' });
const main = h('main');

function show(id: string) {
  if (!tabs.some((t) => t.id === id)) id = 'qpcr';
  for (const t of tabs) {
    if (t.id === id && !panels.has(id)) {
      const el = h('div', { role: 'tabpanel' });
      t.mount(el);
      panels.set(id, el);
      main.append(el);
    }
    panels.get(t.id)?.toggleAttribute('hidden', t.id !== id);
  }
  nav.querySelectorAll('button').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === id)));
  if (location.hash !== `#${id}`) history.replaceState(null, '', `#${id}`);
}

for (const t of tabs) nav.append(h('button', { role: 'tab', 'data-tab': t.id, onclick: () => show(t.id) }, t.label));

document.getElementById('app')!.append(
  h('header', null, h('div', { class: 'brand' }, h('h1', null, 'labcalc'), h('span', null, 'qPCR and cell viability calculators')), h('div', { class: 'privacy' }, '🔒 Runs entirely in your browser. No data is uploaded.')),
  nav,
  main,
  h('footer', null, 'labcalc · free to use · calculations validated against SciPy · see Methods for formulas and references'),
);
window.addEventListener('hashchange', () => show(location.hash.slice(1)));
show(location.hash.slice(1));
