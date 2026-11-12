#!/usr/bin/env python3
"""
fig96_export_tiff.py - the production raster for every figure, beside the vector.

    python3 src/fig96_export_tiff.py

Run it after the figures build. For each figure in results/figures it writes
<stem>.tif beside the existing <stem>.pdf and <stem>.png, then checks the
result against what a production department can accept.

WHY TIFF AND WHY FLATTENED. The figure scripts write a vector PDF and a 600 dpi
PNG. The PDF is what production should set from, and it is supplied. The PNG is
RGBA: matplotlib writes an alpha channel whether or not anything is
transparent, and an alpha channel in a submitted raster is resolved by whatever
opens the file, so a panel drawn on transparency can arrive on black. Each TIFF
is therefore composited onto opaque white and written as 3-channel RGB.

Conversion is pixel-exact: the PNG is lossless at 600 dpi, so this step
re-containers the same pixels rather than re-rendering, and the physical size
is carried over from the figure's own dpi tag. Nothing is resampled.

CHECKS, all of which must pass:
  - the TIFF carries the same pixel dimensions as the PNG it came from
  - mode is RGB, with no alpha channel left
  - the resolution tag is present and at least 300 dpi
  - the physical width agrees with the figure's PDF page width to within 0.2 mm
  - the file is smaller than 50 MB, the usual per-file submission ceiling

OUTPUT
  results/figures/<stem>.tif                one per figure
  logs/fig96_export_tiff.log
"""
from __future__ import annotations

import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)

from atlas_common import DIR_FIGS, Tee, ensure_dirs  # noqa: E402

from PIL import Image  # noqa: E402

MM_PER_IN = 25.4
MIN_DPI = 300
MAX_BYTES = 50 * 1024 * 1024
TOL_MM = 0.2

FAILED = []


def pdf_page_mm(path):
    """The page box of a one-page PDF, in millimetres, read from the file.

    No PDF library: the figure PDFs are written by matplotlib and carry a
    single /MediaBox in clear text, so the page size can be read without
    adding a dependency to a step whose only job is a file conversion.
    """
    with open(path, "rb") as fh:
        head = fh.read(4096).decode("latin-1", "replace")
    m = re.search(r"/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+"
                  r"([\d.+-]+)\s+([\d.+-]+)\s*\]", head)
    if not m:
        return None, None
    x0, y0, x1, y1 = (float(g) for g in m.groups())
    return ((x1 - x0) * MM_PER_IN / 72.0, (y1 - y0) * MM_PER_IN / 72.0)


def main():
    ensure_dirs()
    log = Tee("fig96_export_tiff")
    log("=== fig96  export a production TIFF for every figure ===\n")

    pngs = sorted(p for p in glob.glob(os.path.join(DIR_FIGS, "*.png"))
                  if not os.path.basename(p).startswith("._"))
    if not pngs:
        sys.exit("FATAL: no PNG in %s. Build the figures first with\n"
                 "    python3 src/fig99_build_and_verify.py" % DIR_FIGS)

    log("%-42s %11s %7s %9s %9s" % ("figure", "pixels", "dpi", "width mm",
                                    "size"))
    for png in pngs:
        stem = os.path.splitext(os.path.basename(png))[0]
        pdf = os.path.join(DIR_FIGS, stem + ".pdf")
        tif = os.path.join(DIR_FIGS, stem + ".tif")

        src = Image.open(png)
        dpi = src.info.get("dpi")
        if not dpi:
            FAILED.append("%s: the PNG carries no resolution tag" % stem)
            continue
        dpi_x, dpi_y = (round(float(d), 1) for d in dpi)

        # Composite onto opaque white rather than dropping the alpha channel,
        # so a partly transparent patch keeps the colour it was drawn with.
        if src.mode in ("RGBA", "LA") or "transparency" in src.info:
            flat = Image.new("RGB", src.size, (255, 255, 255))
            flat.paste(src.convert("RGBA"), mask=src.convert("RGBA").split()[-1])
        else:
            flat = src.convert("RGB")

        flat.save(tif, format="TIFF", compression="tiff_lzw",
                  dpi=(dpi_x, dpi_y))

        out = Image.open(tif)
        w_mm = out.width / dpi_x * MM_PER_IN
        h_mm = out.height / dpi_y * MM_PER_IN
        nbytes = os.path.getsize(tif)
        log("%-42s %5dx%-5d %7.0f %9.1f %7.1f MB"
            % (stem + ".tif", out.width, out.height, dpi_x, w_mm,
               nbytes / 1e6))

        if out.size != src.size:
            FAILED.append("%s: %s pixels in, %s out" % (stem, src.size,
                                                        out.size))
        if out.mode != "RGB":
            FAILED.append("%s: mode is %s, not RGB" % (stem, out.mode))
        if min(dpi_x, dpi_y) < MIN_DPI:
            FAILED.append("%s: %.0f dpi, below the %d dpi floor"
                          % (stem, min(dpi_x, dpi_y), MIN_DPI))
        if nbytes > MAX_BYTES:
            FAILED.append("%s: %.1f MB, over the %d MB ceiling"
                          % (stem, nbytes / 1e6, MAX_BYTES // (1024 * 1024)))
        if os.path.exists(pdf):
            pw, ph = pdf_page_mm(pdf)
            if pw is None:
                FAILED.append("%s: no /MediaBox in the PDF" % stem)
            else:
                if abs(pw - w_mm) > TOL_MM or abs(ph - h_mm) > TOL_MM:
                    FAILED.append(
                        "%s: TIFF is %.2f x %.2f mm, PDF is %.2f x %.2f mm"
                        % (stem, w_mm, h_mm, pw, ph))
        else:
            FAILED.append("%s: no vector PDF beside it" % stem)

    log("\n%d TIFF(s) written to %s" % (len(pngs), DIR_FIGS))
    if FAILED:
        sys.exit("FATAL: %d check(s) failed:\n  %s"
                 % (len(FAILED), "\n  ".join(FAILED)))
    log("every TIFF matches its vector PDF in size and is RGB at %d dpi or "
        "better." % MIN_DPI)


if __name__ == "__main__":
    main()
