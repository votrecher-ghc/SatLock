"""Crop PDF whitespace while preserving equal page sizes within a figure group.

The MATLAB pipeline exports vector PDFs first. This script measures the visible
content in each one-page PDF, then rewrites all PDFs in the group to a shared
tight page size. It keeps vector content by embedding the original page with a
clip rectangle instead of rasterizing the final output.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import fitz  # PyMuPDF
import numpy as np


def visible_bbox(page: fitz.Page, dpi: int = 216, threshold: int = 248) -> fitz.Rect:
    scale = dpi / 72.0
    pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, pix.n)
    rgb = arr[:, :, : min(3, pix.n)]
    mask = np.any(rgb < threshold, axis=2)
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return fitz.Rect(page.rect)
    return fitz.Rect(
        xs.min() / scale,
        ys.min() / scale,
        (xs.max() + 1) / scale,
        (ys.max() + 1) / scale,
    )


def padded_clip(page_rect: fitz.Rect, bbox: fitz.Rect, pad: float) -> fitz.Rect:
    """Return a safe clip that includes the full visible bbox plus padding."""
    return fitz.Rect(
        max(page_rect.x0, bbox.x0 - pad),
        max(page_rect.y0, bbox.y0 - pad),
        min(page_rect.x1, bbox.x1 + pad),
        min(page_rect.y1, bbox.y1 + pad),
    )


def centered_destination(page_width: float, page_height: float, clip: fitz.Rect) -> fitz.Rect:
    x0 = (page_width - clip.width) / 2.0
    y0 = (page_height - clip.height) / 2.0
    return fitz.Rect(x0, y0, x0 + clip.width, y0 + clip.height)


def normalize_group(paths: list[Path], pad: float) -> None:
    page_rects: list[fitz.Rect] = []
    bboxes: list[fitz.Rect] = []
    for path in paths:
        with fitz.open(path) as doc:
            page = doc[0]
            page_rects.append(fitz.Rect(page.rect))
            bboxes.append(visible_bbox(page))

    clips = [padded_clip(page_rect, bbox, pad) for page_rect, bbox in zip(page_rects, bboxes)]
    target_w = max(clip.width for clip in clips) + 2.0 * pad
    target_h = max(clip.height for clip in clips) + 2.0 * pad

    for path, clip in zip(paths, clips):
        src_doc = fitz.open(path)
        out_doc = fitz.open()
        out_page = out_doc.new_page(width=target_w, height=target_h)
        out_page.show_pdf_page(centered_destination(target_w, target_h, clip), src_doc, 0, clip=clip)

        tmp_path = path.with_suffix(path.suffix + ".tmp")
        out_doc.save(tmp_path, garbage=4, deflate=True)
        out_doc.close()
        src_doc.close()
        os.replace(tmp_path, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pad", type=float, default=3.0)
    parser.add_argument("pdfs", nargs="+", type=Path)
    args = parser.parse_args()
    normalize_group(args.pdfs, args.pad)


if __name__ == "__main__":
    main()
