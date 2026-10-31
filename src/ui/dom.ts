import { decodeBytes } from '../lib/parse';

// Tiny DOM helper: h('div', { class: 'x', onclick: fn }, child1, child2)
type Child = Node | string | number | null | undefined | false;
type Attrs = Record<string, unknown>;

export function h<K extends keyof HTMLElementTagNameMap>(tag: K, attrs: Attrs | null = null, ...children: (Child | Child[])[]): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  if (attrs)
    for (const [k, v] of Object.entries(attrs)) {
      if (v === undefined || v === null || v === false) continue;
      if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v as EventListener);
      else if (k === 'value' || k === 'checked' || k === 'selected' || k === 'disabled') (el as unknown as Record<string, unknown>)[k] = v;
      else if (k === 'html') el.innerHTML = String(v);
      else el.setAttribute(k, v === true ? '' : String(v));
    }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    el.append(c instanceof Node ? c : String(c));
  }
  return el;
}

export function fmt(v: number, digits = 3): string {
  if (!Number.isFinite(v)) return '';
  if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e5)) return v.toExponential(2);
  return v.toFixed(digits);
}

export function table(header: string[], rows: Child[][], cls = ''): HTMLTableElement {
  return h('table', { class: `data ${cls}` }, h('thead', null, h('tr', null, header.map((x) => h('th', null, x)))), h('tbody', null, rows.map((r) => h('tr', null, r.map((c) => h('td', null, c))))));
}

export function card(title: string, ...body: (Child | Child[])[]): HTMLElement {
  return h('section', { class: 'card' }, h('h2', null, title), ...body);
}

/** Read a File as text, handling UTF-16 exports. */
export async function readFileText(file: File): Promise<string> {
  return decodeBytes(await file.arrayBuffer());
}

/**
 * Re-render a container while keeping keyboard focus and scroll positions:
 * inputs carry a data-k key so the focused field can be found again.
 */
export function rerender(container: HTMLElement, build: () => Node): void {
  const active = document.activeElement as HTMLElement | null;
  const key = active?.dataset?.k;
  const scrolls = [...container.querySelectorAll<HTMLElement>('[data-scroll]')].map((e) => [e.dataset.scroll!, e.scrollTop] as const);
  const y = window.scrollY;
  container.replaceChildren(build());
  for (const [k, top] of scrolls) {
    const e = container.querySelector<HTMLElement>(`[data-scroll="${k}"]`);
    if (e) e.scrollTop = top;
  }
  if (key) container.querySelector<HTMLElement>(`[data-k="${CSS.escape(key)}"]`)?.focus();
  window.scrollTo(0, y);
}
