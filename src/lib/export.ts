// Client-side downloads: CSV tables, SVG and PNG figures. Nothing is uploaded.

export type Cell = string | number | null | undefined;

function csvCell(v: Cell): string {
  if (v === null || v === undefined) return '';
  if (typeof v === 'number') return Number.isFinite(v) ? String(v) : '';
  return /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;
}

export function toCsv(header: string[], rows: Cell[][]): string {
  return [header, ...rows].map((r) => r.map(csvCell).join(',')).join('\r\n') + '\r\n';
}

export function download(filename: string, data: Blob | string, type = 'text/plain'): void {
  const blob = typeof data === 'string' ? new Blob([data], { type: `${type};charset=utf-8` }) : data;
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function downloadSvg(filename: string, svg: string): void {
  download(filename, svg, 'image/svg+xml');
}

/** Rasterise an SVG string to PNG at the given scale (3 = about 300 dpi for a 100 dpi layout). */
export async function downloadPng(filename: string, svg: string, scale = 3): Promise<void> {
  const w = Number(svg.match(/width="(\d+(?:\.\d+)?)"/)?.[1] ?? 800);
  const h = Number(svg.match(/height="(\d+(?:\.\d+)?)"/)?.[1] ?? 600);
  const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }));
  try {
    const img = new Image();
    await new Promise<void>((res, rej) => {
      img.onload = () => res();
      img.onerror = () => rej(new Error('Could not render the figure.'));
      img.src = url;
    });
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(w * scale);
    canvas.height = Math.round(h * scale);
    const ctx = canvas.getContext('2d')!;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0, w, h);
    const blob = await new Promise<Blob | null>((res) => canvas.toBlob(res, 'image/png'));
    if (blob) download(filename, blob);
  } finally {
    URL.revokeObjectURL(url);
  }
}
