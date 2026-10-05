#!/usr/bin/env python3
"""Draw the graphical abstract for the geoExposeR SoftwareX paper.

Three panels follow the package workflow: combine modeled and measured exposure
probabilities with population weights, impute each unit's exposure and fit a
regression model to every imputed dataset, and pool the exposure effects with
Rubin's rules. All shapes are drawn directly in SVG (no third-party icons).

Elsevier asks for at least 1328 x 531 px (w x h), readable at 13 x 5 cm, so the
canvas is 1328 units wide and every label is at least 26 units tall (about 7 pt
at 13 cm). Exposure categories use one blue ramp, light to dark for low to high,
and are named in the legend, so color is never the only cue.

Outputs:
  man/figures/graphical_abstract.png       2656 x 1062 px, for the README
  paper/geoExposeR_graphical_abstract.svg  editable vector source
  paper/geoExposeR_graphical_abstract.pdf  vector, fonts embedded, 13 x 5.2 cm
  paper/geoExposeR_graphical_abstract.tif  5312 x 2124 px (4x minimum), LZW compressed
  paper/geoExposeR_graphical_abstract.png  same size as the TIFF
The paper/ manuscript folder is not part of the repository, so the submission
files are written only when a local paper/ folder exists.

Requirements: Python 3 with Pillow, and R with the rsvg package.
Usage (from the package root): python3 scripts/generate_graphical_abstract.py
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image, ImageFont

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper"
NAME = "geoExposeR_graphical_abstract"
README_PNG = ROOT / "man" / "figures" / "graphical_abstract.png"
W, H = 1328, 531    # Elsevier minimum size in px (w x h)
WIDTH_CM = 13.0     # printed width at which the labels must stay readable
SCALE = 4           # the TIFF is SCALE times the minimum size

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
HAIRLINE, PANEL, PANEL_EDGE = "#c3c2b7", "#f6f5f2", "#e1e0d9"
# Exposure categories low, middle and high: ordinal steps of one blue ramp
RAMP = ("#86b6ef", "#2a78d6", "#104281")
FONT = "Arial, Helvetica, sans-serif"
FONT_FILES = {
    "normal": "/System/Library/Fonts/Supplemental/Arial.ttf",
    "bold": "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "italic": "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
}

# Panel geometry: three panels with arrow gaps between them; the first is wider
# so that its three bar labels do not touch
PY, PH, GAP, PAD = 18, 440, 56, 24
PWS = (420, 380, 380)
PX = (18, 18 + PWS[0] + GAP, 18 + PWS[0] + PWS[1] + 2 * GAP)
SPACE_EM = 0.278  # Arial space advance; spaces are written as tspan offsets

svg = []
texts = []  # (string, size, face, max_width) for the fit check


def n(v):
    """Format a coordinate compactly."""
    return f"{v:.2f}".rstrip("0").rstrip(".")


def face(weight, style):
    return "bold" if weight == "bold" else ("italic" if style == "italic" else "normal")


def text_width(s, size, key="normal"):
    path = FONT_FILES[key]
    if Path(path).exists():
        return ImageFont.truetype(path, size).getlength(s)
    return 0.55 * size * len(s)


def text(x, y, content, size=26, anchor="start", weight="normal", style="normal",
         fill=INK, max_width=None):
    """Add a text element; content is a string or a list of (string, attrs)."""
    attrs = f'x="{n(x)}" y="{n(y)}" font-size="{size}" fill="{fill}"'
    if anchor != "start":
        attrs += f' text-anchor="{anchor}"'
    if weight != "normal":
        attrs += f' font-weight="{weight}"'
    if style != "normal":
        attrs += f' font-style="{style}"'
    if isinstance(content, str):
        content = [(content, {})]
    # Each word is its own tspan and each space becomes a dx offset: cairo
    # embeds the space glyph as a Type 3 font, which PDF preflight checks flag
    body, pending = "", 0.0
    for s, extra in content:
        seg_size = float(extra.get("font-size", size))
        for word in s.split(" "):
            if word:
                attrs_w = dict(extra)
                if pending:
                    attrs_w["dx"] = n(pending)
                ex = " ".join(f'{k}="{v}"' for k, v in attrs_w.items())
                body += f"<tspan {ex}>{escape(word)}</tspan>" if ex else f"<tspan>{escape(word)}</tspan>"
                pending = 0.0
            pending += SPACE_EM * seg_size
        pending -= SPACE_EM * seg_size  # split() adds one boundary per segment
    svg.append(f"<text {attrs}>{body}</text>")
    texts.append(("".join(s for s, _ in content), size, face(weight, style), max_width))


def arrow(x1, y1, x2, y2, color=INK2, width=4, head=13):
    """Straight arrow from (x1, y1) to (x2, y2) with a filled head."""
    dx, dy = x2 - x1, y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / length, dy / length
    bx, by = x2 - ux * head, y2 - uy * head
    svg.append(
        f'<line x1="{n(x1)}" y1="{n(y1)}" x2="{n(bx)}" y2="{n(by)}" '
        f'stroke="{color}" stroke-width="{width}" stroke-linecap="round"/>'
    )
    px, py = -uy * head * 0.6, ux * head * 0.6
    svg.append(
        f'<polygon points="{n(x2)},{n(y2)} {n(bx + px)},{n(by + py)} '
        f'{n(bx - px)},{n(by - py)}" fill="{color}"/>'
    )


def stacked_bar(cx, top, base, width, probs, gap=2, radius=4):
    """Vertical stacked bar: lowest category on the baseline, 2-unit gaps."""
    x0 = cx - width / 2
    usable = (base - top) - gap * (len(probs) - 1)
    y = base
    for i, (p, color) in enumerate(zip(probs, RAMP)):
        h = usable * p
        y_top = y - h
        if i == len(probs) - 1:
            # Round only the free end; the baseline end stays square
            r = min(radius, h / 2)
            svg.append(
                f'<path d="M{n(x0)},{n(y)} V{n(y_top + r)} Q{n(x0)},{n(y_top)} '
                f'{n(x0 + r)},{n(y_top)} H{n(x0 + width - r)} Q{n(x0 + width)},'
                f'{n(y_top)} {n(x0 + width)},{n(y_top + r)} V{n(y)} Z" fill="{color}"/>'
            )
        else:
            svg.append(
                f'<rect x="{n(x0)}" y="{n(y_top)}" width="{n(width)}" '
                f'height="{n(h)}" fill="{color}"/>'
            )
        y = y_top - gap


def tile(x, y, cats, w=128, h=98, pad=7, gap=3):
    """A map of geographic units, each colored by its imputed category."""
    svg.append(
        f'<rect x="{n(x)}" y="{n(y)}" width="{w}" height="{h}" rx="8" '
        f'fill="#ffffff" stroke="{HAIRLINE}" stroke-width="1.5"/>'
    )
    rows, cols = len(cats), len(cats[0])
    cw = (w - 2 * pad - gap * (cols - 1)) / cols
    ch = (h - 2 * pad - gap * (rows - 1)) / rows
    for r, row in enumerate(cats):
        for c, k in enumerate(row):
            svg.append(
                f'<rect x="{n(x + pad + c * (cw + gap))}" y="{n(y + pad + r * (ch + gap))}" '
                f'width="{n(cw)}" height="{n(ch)}" rx="3" fill="{RAMP[k]}"/>'
            )


def person(cx, top, fill=INK2):
    svg.append(f'<circle cx="{n(cx)}" cy="{n(top + 6)}" r="6" fill="{fill}"/>')
    svg.append(
        f'<rect x="{n(cx - 8)}" y="{n(top + 14)}" width="16" height="14" rx="6" fill="{fill}"/>'
    )


def interval(y, lo, hi, est, color, width, r):
    svg.append(
        f'<line x1="{n(lo)}" y1="{n(y)}" x2="{n(hi)}" y2="{n(y)}" stroke="{color}" '
        f'stroke-width="{width}" stroke-linecap="round"/>'
    )
    svg.append(
        f'<circle cx="{n(est)}" cy="{n(y)}" r="{r}" fill="{color}" '
        f'stroke="#ffffff" stroke-width="2"/>'
    )


def panel(x, pw, number, title):
    svg.append(
        f'<rect x="{x}" y="{PY}" width="{pw}" height="{PH}" rx="14" fill="{PANEL}" '
        f'stroke="{PANEL_EDGE}" stroke-width="1.5"/>'
    )
    svg.append(f'<circle cx="{x + 38}" cy="{PY + 40}" r="19" fill="{INK}"/>')
    text(x + 38, PY + 49, str(number), size=26, anchor="middle", weight="bold",
         fill="#ffffff")
    text(x + 68, PY + 50, title, size=31, weight="bold", max_width=pw - 68 - PAD)


def panel_combine(x, pw):
    panel(x, pw, 1, "Combine sources")
    inner = pw - 2 * PAD
    mid = x + pw / 2
    step = 130
    centers = (mid - step, mid, mid + step)
    w = 0.4  # illustrative private-well share
    modeled, measured = (0.55, 0.30, 0.15), (0.85, 0.12, 0.03)
    combined = tuple(w * a + (1 - w) * b for a, b in zip(modeled, measured))
    for cx, probs in zip(centers, (modeled, measured, combined)):
        stacked_bar(cx, 104, 284, 56, probs)
    for cx, op in zip(((centers[0] + centers[1]) / 2, (centers[1] + centers[2]) / 2), "+="):
        text(cx, 208, op, size=40, anchor="middle", fill=INK2)
    for cx, label in zip(centers, ("Modeled", "Measured", "Combined")):
        text(cx, 320, label, anchor="middle", max_width=step - 6)
    text(centers[0], 352, "× w", anchor="middle", style="italic", fill=INK2)
    text(centers[1], 352, "× (1 − w)", anchor="middle", style="italic", fill=INK2)
    text(mid, 404, [("w", {"font-style": "italic"}), (" = population share", {})],
         anchor="middle", fill=INK2, max_width=inner)
    text(mid, 434, "served by the modeled source", anchor="middle", fill=INK2,
         max_width=inner)


def panel_impute(x, pw):
    panel(x, pw, 2, "Impute and fit")
    left = x + PAD
    draws = (
        ((0, 0, 1, 0), (0, 2, 0, 0), (1, 0, 0, 0)),
        ((0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 2, 0)),
        ((0, 0, 0, 1), (1, 0, 0, 0), (0, 2, 0, 0)),
    )
    for k, cats in enumerate(draws):
        tile(left + 15 * k, 96 + 15 * k, cats)
    tiles_right = left + 30 + 128
    text(tiles_right + 22, 162, [("m", {"font-style": "italic"}), (" imputed", {})],
         size=27, max_width=x + pw - PAD - tiles_right - 22)
    text(tiles_right + 22, 196, "datasets", size=27)
    ax = left + 72  # under the tile stack
    arrow(ax, 234, ax, 280)
    for k in range(3):
        person(ax + 28 + 18 * k, 238)
    label_x = ax + 84
    text(label_x, 264, "health records", fill=INK2, max_width=x + pw - PAD - label_x)
    box_x, box_w = left, pw - 2 * PAD
    svg.append(
        f'<rect x="{box_x}" y="284" width="{box_w}" height="134" rx="10" '
        f'fill="#ffffff" stroke="{HAIRLINE}" stroke-width="1.5"/>'
    )
    cx = box_x + box_w / 2
    text(cx, 322, "Regression per dataset", size=27, anchor="middle", weight="bold",
         max_width=box_w - 24)
    text(cx, 358, "fixed, random, or mixed", anchor="middle", fill=INK2, max_width=box_w - 24)
    text(cx, 392, "linear, GLM, multinomial", anchor="middle", fill=INK2, max_width=box_w - 24)


def panel_pool(x, pw):
    panel(x, pw, 3, "Pool estimates")
    left = x + PAD
    inner = pw - 2 * PAD
    text(left, 112, "Effect vs low exposure", fill=INK2, max_width=inner)
    zero = x + 186
    svg.append(
        f'<line x1="{zero}" y1="124" x2="{zero}" y2="250" stroke="{HAIRLINE}" '
        f'stroke-width="2.5" stroke-dasharray="6 6"/>'
    )
    text(zero, 276, "0", anchor="middle", fill=INK2)
    rows = (("Middle", 152, (-38, 86, 24), (-64, 114, 26)),
            ("High", 218, (12, 126, 68), (-12, 154, 71)))
    for label, cy, one, pooled in rows:
        text(left, cy + 9, label, size=27)
        interval(cy - 11, zero + one[0], zero + one[1], zero + one[2], MUTED, 3, 6)
        interval(cy + 11, zero + pooled[0], zero + pooled[1], zero + pooled[2], INK, 4, 7)
    text(zero + 56, 304, "Exposure effect", anchor="middle", fill=INK2)
    for y, color, width, r, label in ((334, MUTED, 3, 6, "one imputed dataset"),
                                      (364, INK, 4, 7, "pooled (Rubin’s rules)")):
        interval(y, left + 2, left + 40, left + 21, color, width, r)
        text(left + 54, y + 8, label, max_width=inner - 54)
    text(left, 412, "Pooled intervals include", weight="bold", max_width=inner)
    text(left, 441, "exposure uncertainty", weight="bold", max_width=inner)


def footer():
    y = 503
    x = PX[0] + 6
    text(x, y, "Exposure", fill=INK2)
    x += text_width("Exposure", 26) + 22
    for color, label in zip(RAMP, ("Low", "Middle", "High")):
        svg.append(f'<rect x="{n(x)}" y="{y - 21}" width="24" height="24" rx="4" fill="{color}"/>')
        text(x + 32, y, label, fill=INK2)
        x += 32 + text_width(label, 26) + 24
    text(PX[2] + PWS[2], y + 1, [("geoExposeR", {"font-weight": "bold", "fill": INK, "font-size": "30"}),
                            ("  ·  arsenic, nitrate, uranium, PM2.5", {})],
         size=28, anchor="end", fill=INK2)


def build():
    svg.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
    panel_combine(PX[0], PWS[0])
    panel_impute(PX[1], PWS[1])
    panel_pool(PX[2], PWS[2])
    mid_y = PY + PH / 2
    for x in (PX[0] + PWS[0], PX[1] + PWS[1]):
        arrow(x + 10, mid_y, x + GAP - 10, mid_y)
    footer()
    height_cm = WIDTH_CM * H / W
    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_CM}cm" '
        f'height="{height_cm:.3f}cm" viewBox="0 0 {W} {H}" font-family="{FONT}">'
        "<title>geoExposeR graphical abstract</title>"
        "<desc>Modeled and measured exposure probabilities are combined with "
        "population weights, each unit's exposure category is imputed many times "
        "and a regression model is fit to each imputed dataset, and the exposure "
        "effects are pooled with Rubin's rules so that their intervals include "
        "exposure uncertainty. Worked examples cover arsenic, nitrate, uranium, "
        "and fine particulate matter (PM2.5).</desc>"
    )
    return head + "".join(svg) + "</svg>\n"


def check_fit():
    ok = True
    for s, size, key, max_width in texts:
        if max_width is not None and text_width(s, size, key) > max_width:
            ok = False
            print(f"warning: '{s}' is {text_width(s, size, key):.0f} units wide, "
                  f"more than its {max_width:.0f}", file=sys.stderr)
    return ok


def main():
    svg_text = build()
    check_fit()
    with tempfile.TemporaryDirectory() as tmp:
        svg_file, pdf_file, png_file = (Path(tmp) / f"{NAME}{ext}"
                                        for ext in (".svg", ".pdf", ".png"))
        svg_file.write_text(svg_text, encoding="utf-8")
        # PDF page in points so that it prints WIDTH_CM wide
        pt_w = WIDTH_CM / 2.54 * 72
        script = (
            f'rsvg::rsvg_pdf("{svg_file}", "{pdf_file}", '
            f'width = {pt_w:.2f}, height = {pt_w * H / W:.2f}); '
            f'rsvg::rsvg_png("{svg_file}", "{png_file}", width = {W * SCALE})'
        )
        subprocess.run(["Rscript", "-e", script], check=True)
        img = Image.open(png_file).convert("RGB")
        img.resize((2 * W, 2 * H), Image.LANCZOS).save(README_PNG, optimize=True)
        print(f"Wrote {README_PNG.relative_to(ROOT)}")
        if not PAPER.is_dir():
            print("No paper/ folder: skipped the submission files")
            return
        dpi = W * SCALE / (WIDTH_CM / 2.54)
        img.save(PAPER / f"{NAME}.tif", compression="tiff_lzw", dpi=(dpi, dpi))
        img.save(PAPER / f"{NAME}.png", optimize=True, dpi=(dpi, dpi))
        shutil.copy(svg_file, PAPER / f"{NAME}.svg")
        shutil.copy(pdf_file, PAPER / f"{NAME}.pdf")
        print(f"Wrote paper/{NAME}.svg, .pdf, .tif and .png "
              f"({img.width} x {img.height} px)")


if __name__ == "__main__":
    main()
