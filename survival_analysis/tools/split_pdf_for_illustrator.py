"""Split R multi-page figure PDFs into one file per page for Illustrator editing.

Why: Illustrator's "PDF Import Options" dialog appears only for multi-page PDFs (or
linked ones). For those, Illustrator creates a new document and imports the
pages as LINKED objects - one non-editable object per page ("a whole image",
no selectable parts, no editable text). A single-page PDF is imported directly
as embedded artwork: every path and every text frame stays editable, no dialog.

Every delivered figure keeps its multi-page PDF; this tool additionally writes
one single-page PDF per page into <figure>/<EXPR>/editable_1page/.

Usage:
    python split_pdf_for_illustrator.py <project_root>
"""
import os
import re
import shutil
import sys

import fitz  # PyMuPDF: page labels / blank-page detection
import pikepdf


def page_label(page):
    """Short, human-readable label for one figure page."""
    words = page.get_text('words')
    cand = [w for w in words if not set(w[4]) <= set('+- ')]
    cand = sorted(cand, key=lambda w: w[1])[:3]
    toks = [w[4] for w in cand]
    if not toks:
        return None

    cutoff = None
    for t in toks:
        m = re.match(r'^\(([A-Za-z0-9_]+)', t)
        if m:
            cutoff = m.group(1)
            break

    first = toks[0]
    if first in ('FPKM', 'TPM'):                       # C-index comparison pages
        return cutoff or first
    if first.endswith(':') and len(toks) > 1:          # "Changhai: HSPD1" + "(Median)"
        parts = [first[:-1], re.sub(r'[^A-Za-z0-9_]', '', toks[1])]
        if cutoff:
            parts.append(cutoff)
        return '_'.join(p for p in parts if p)
    return re.sub(r'[^A-Za-z0-9_.+-]', '', first)[:48] or None


def split(pdf_path, out_dir):
    src = fitz.open(pdf_path)
    if src.page_count < 2:
        src.close()
        return []
    stem = os.path.splitext(os.path.basename(pdf_path))[0]
    os.makedirs(out_dir, exist_ok=True)
    pdf = pikepdf.open(pdf_path)
    written, idx = [], 0
    for i in range(src.page_count):
        if not src[i].get_text().strip() and len(src[i].get_drawings()) <= 2:
            continue                                   # blank filler page
        idx += 1
        label = page_label(src[i]) or ('page%02d' % (i + 1))
        name = '%s_%02d_%s.pdf' % (stem, idx, label)
        dst = pikepdf.new()
        dst.pages.append(pdf.pages[i])
        dst.save(os.path.join(out_dir, name),
                 compress_streams=True,
                 object_stream_mode=pikepdf.ObjectStreamMode.disable)
        dst.close()
        written.append(name)
    pdf.close()
    src.close()
    return written


if __name__ == '__main__':
    root = sys.argv[1]
    total = 0
    for expr in ('FPKM', 'TPM'):
        fig_dir = os.path.join(root, 'figure', expr)
        if not os.path.isdir(fig_dir):
            continue
        out_dir = os.path.join(fig_dir, 'editable_1page')
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)
        for f in sorted(os.listdir(fig_dir)):
            if not f.lower().endswith('.pdf'):
                continue
            got = split(os.path.join(fig_dir, f), out_dir)
            total += len(got)
            if got:
                print('%-56s -> %2d single-page files' % (f, len(got)))
                for g in got[:3]:
                    print('        ', g)
                if len(got) > 3:
                    print('         ... +%d more' % (len(got) - 3))
    print('\ntotal single-page PDFs written:', total)
