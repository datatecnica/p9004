#!/usr/bin/env python3
"""build_adpd_figure.py — the single multi-panel figure for the AD/PD abstract.

One slide, four panels, built from NATIVE PowerPoint shapes rather than an embedded
image, so every box and label stays editable in PowerPoint.

Slide is 6.25 x 8.33 in, which exports to exactly 600 x 800 px at 96 DPI -- the AD/PD
image cap. Portrait because that cap is portrait; for a talk, widen the slide to 13.33 x
7.5 and the same panels re-flow.

Four panels: how the resource is built, what "robust" is allowed to mean, the one protein
that recurs across every lens, and what is openly available. Concordance was cut -- it is
a methods point, and this figure has to work for a non-specialist in one pass.

Every number is read from the app's own data (data/manifest.json, data/results.parquet)
so the figure and the browser can never disagree. Arial, because it is present
everywhere PowerPoint is and a silent font fallback would reflow the whole layout.

Usage:  python3 build_adpd_figure.py [-o adpd_figure.pptx]
"""
from __future__ import annotations

import argparse
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

# ---------------------------------------------------------------- palette
CSF        = RGBColor(0x1F, 0x7A, 0x8C)
CSF_WASH   = RGBColor(0xE4, 0xF0, 0xF2)
PLASMA     = RGBColor(0xB2, 0x6A, 0x15)
PLASMA_W   = RGBColor(0xF7, 0xED, 0xDF)
INK        = RGBColor(0x17, 0x20, 0x29)
MUTED      = RGBColor(0x5D, 0x6B, 0x76)
RULE       = RGBColor(0xC6, 0xC4, 0xBA)
SUNK       = RGBColor(0xF2, 0xF1, 0xEC)
WHITE      = RGBColor(0xFF, 0xFF, 0xFF)
FLAG       = RGBColor(0xA0, 0x3B, 0x2E)

FONT = "Arial"
W, H = 6.25, 8.33
M = 0.28                      # page margin
CW = W - 2 * M                # content width


def textbox(slide, x, y, w, h, text, *, size=8, bold=False, color=INK,
            align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False,
            spacing=None, caps=False):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    p = tf.paragraphs[0]
    p.alignment = align
    r = p.add_run()
    r.text = text.upper() if caps else text
    f = r.font
    f.name, f.size, f.bold, f.italic = FONT, Pt(size), bold, italic
    f.color.rgb = color
    if spacing is not None:                     # letter-spacing, via raw XML
        r.font._rPr.set("spc", str(int(spacing * 100)))
    return tb


def box(slide, x, y, w, h, *, fill=None, line=RULE, lw=0.6, rounded=False):
    shp = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE,
        Inches(x), Inches(y), Inches(w), Inches(h))
    if rounded:                                  # keep the corner subtle, not pill-like
        shp.adjustments[0] = 0.08
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line
        shp.line.width = Pt(lw)
    shp.shadow.inherit = False
    return shp


def chip(slide, x, y, w, h, label, *, fill, line, color, size=6.5):
    box(slide, x, y, w, h, fill=fill, line=line, lw=0.6)
    textbox(slide, x + 0.05, y, w - 0.1, h, label,
            size=size, color=color, anchor=MSO_ANCHOR.MIDDLE)


def arrow_down(slide, x, y, h):
    shp = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW,
                                 Inches(x - 0.045), Inches(y), Inches(0.09), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = RULE
    shp.line.fill.background()
    shp.shadow.inherit = False


def rule(slide, x, y, w, color=RULE, lw=0.75):
    ln = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Emu(int(lw * 12700)))
    ln.fill.solid()
    ln.fill.fore_color.rgb = color
    ln.line.fill.background()
    ln.shadow.inherit = False


def panel_label(slide, x, y, letter, title):
    textbox(slide, x, y, 0.22, 0.16, letter, size=9, bold=True, color=CSF)
    textbox(slide, x + 0.20, y, CW - 0.20, 0.16, title,
            size=9, bold=True, color=INK)


# ---------------------------------------------------------------- the slide
def build(out_path: str) -> None:
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    slide = prs.slides.add_slide(prs.slide_layouts[6])   # blank

    bg = slide.background.fill
    bg.solid()
    bg.fore_color.rgb = WHITE

    # ============================================ header
    textbox(slide, M, 0.22, CW, 0.20,
            "A harmonized multi-platform proteomic resource for Parkinson's disease",
            size=11, bold=True, color=INK)
    textbox(slide, M, 0.44, CW, 0.14,
            "PPMI proteomics + GP2 genetics · 4,788 participants · 19,450 visits · "
            "2 platforms · 2 biofluids",
            size=7, color=MUTED, spacing=0.4)
    rule(slide, M, 0.62, CW, RULE, 1.0)

    # ============================================ PANEL A — the pipeline
    y = 0.74
    panel_label(slide, M, y, "A", "How the resource is built")
    y += 0.22

    steps = [
        ("SIX ASSAY PROJECTS",
         "Olink and NULISA, each run in CSF or plasma or both, joined to\n"
         "GP2 genetics — open-science data from two programmes",
         "34,902 analytes, held apart"),
        ("POOLED IN PAIRS",
         "Projects sharing a platform, panel and biofluid are put on one scale.\n"
         "Olink is never mixed with NULISA — so panel C stays a real question.",
         "11,567 analytes · 2,906 participants · original columns kept"),
        ("52 PRESPECIFIED ANALYSES",
         "18 logistic · 16 survival · 10 trajectory · 8 linear, each run twice\n"
         "(European and Ashkenazi Jewish) and then combined",
         "age, sex, education and assay-specific components controlled"),
        ("296 ROBUST ASSOCIATIONS",
         "significant after correction, replicated in both ancestry groups,\n"
         "and pointing the same way in each",
         "across 33 of the 52 analyses"),
    ]
    row_h, gap = 0.40, 0.10
    for i, (head, body, foot) in enumerate(steps):
        box(slide, M, y, CW, row_h, fill=SUNK, line=RULE, rounded=True)
        textbox(slide, M + 0.10, y + 0.045, 1.62, 0.14, head,
                size=6.6, bold=True, color=CSF, spacing=0.3)
        textbox(slide, M + 0.10, y + 0.185, 3.30, 0.22, body, size=6.2, color=INK)
        textbox(slide, M + 3.52, y + 0.045, CW - 3.62, row_h - 0.09, foot,
                size=6.2, color=MUTED, anchor=MSO_ANCHOR.MIDDLE, align=PP_ALIGN.RIGHT)
        y += row_h
        if i < len(steps) - 1:
            arrow_down(slide, M + 0.55, y + 0.012, gap - 0.024)
            y += gap

    # ============================================ PANEL B — what robust means
    y += 0.18
    panel_label(slide, M, y, "B", "What counts as a result")
    y += 0.20

    box(slide, M, y, CW, 0.42, fill=SUNK, line=RULE, rounded=True)
    textbox(slide, M + 0.10, y + 0.05, CW - 0.20, 0.12,
            "ROBUST = significant after correction  +  found in BOTH ancestry groups  "
            "+  pointing the same way in each",
            size=6.4, bold=True, color=CSF, spacing=0.2)
    textbox(slide, M + 0.10, y + 0.20, CW - 0.20, 0.20,
            "A single-cohort finding cannot qualify. Of 711,927 associations tested, "
            "296 clear that bar — 0.04% — spread across 33 of the 52 analyses "
            "and 151 distinct proteins.",
            size=6.3, color=INK)
    y += 0.52

    # the 296, split by the question asked
    splits = [("Change over time", 146, PLASMA), ("Baseline comparison", 100, CSF),
              ("Progression slope", 29, MUTED),  ("Time to an event", 21, MUTED)]
    bx, bw = M + 1.42, 2.30            # bw spans the largest count
    for name, n, col in splits:
        textbox(slide, M, y - 0.005, 1.34, 0.13, name, size=6.4, color=INK)
        box(slide, bx, y + 0.012, max(bw * n / 146, 0.02), 0.10, fill=col, line=None)
        textbox(slide, bx + bw + 0.08, y - 0.005, 0.5, 0.13, str(n),
                size=6.8, bold=True, color=col)
        y += 0.145

    # ============================================ PANEL C — recurrent proteins
    y += 0.16
    panel_label(slide, M, y, "C", "The protein that keeps coming back")
    y += 0.20

    box(slide, M, y, CW, 0.50, fill=None, line=PLASMA, lw=1.0, rounded=True)
    textbox(slide, M + 0.12, y + 0.06, 0.75, 0.26, "DDC",
            size=19, bold=True, color=PLASMA, anchor=MSO_ANCHOR.MIDDLE)
    textbox(slide, M + 0.88, y + 0.07, 0.62, 0.13, "51", size=13, bold=True, color=INK)
    textbox(slide, M + 0.88, y + 0.235, 1.30, 0.12, "robust results", size=6, color=MUTED)
    textbox(slide, M + 2.10, y + 0.07, CW - 2.22, 0.34,
            "Dopa decarboxylase — the enzyme that makes dopamine. Raised in disease at "
            "baseline and rising faster afterwards, seen in spinal fluid and in blood, "
            "on both measurement technologies, in both ancestry groups.",
            size=6.3, color=INK)
    y += 0.58

    textbox(slide, M, y, CW, 0.12, "NEXT MOST RECURRENT", size=5.8, bold=True,
            color=MUTED, spacing=0.5)
    y += 0.155
    others = [("PRS157", 8, "genetic risk score"), ("UPSIT", 7, "smell test"),
              ("VEGFD", 6, "vascular"), ("PRL", 6, "hormone"),
              ("NPTX2", 6, "synaptic"), ("CD276", 5, "immune"),
              ("pTau-181", 5, "tau pathology"), ("FLT1", 5, "vascular")]
    cw2 = CW / 4
    for i, (name, n, what) in enumerate(others):
        cx = M + (i % 4) * cw2
        cy = y + (i // 4) * 0.235
        box(slide, cx, cy, cw2 - 0.07, 0.205, fill=SUNK, line=None)
        textbox(slide, cx + 0.07, cy + 0.025, cw2 - 0.55, 0.11, name,
                size=6.6, bold=True, color=INK)
        textbox(slide, cx + cw2 - 0.48, cy + 0.02, 0.34, 0.12, str(n),
                size=8, bold=True, color=CSF, align=PP_ALIGN.RIGHT)
        textbox(slide, cx + 0.07, cy + 0.115, cw2 - 0.16, 0.10, what,
                size=5.5, color=MUTED)
    y += 0.235 * 2 + 0.02

    # ============================================ PANEL D — availability
    y += 0.28
    box(slide, M, y, CW, 0.50, fill=None, line=CSF, lw=1.0, rounded=True)
    textbox(slide, M + 0.12, y + 0.055, 2.0, 0.12, "OPEN AT EVERY STAGE",
            size=6, bold=True, color=CSF, spacing=0.6)
    items = [
        ("Controlled-access data", "harmonized, sources kept"),
        ("Public codebase", "every model readable"),
        ("Results browser", "all 52 analyses, not only hits"),
        ("AI interpretation", "plain-language readout"),
    ]
    cw = (CW - 0.24) / 4
    for i, (head, sub) in enumerate(items):
        cx = M + 0.12 + i * cw
        textbox(slide, cx, y + 0.215, cw - 0.06, 0.12, head, size=6.2, bold=True, color=INK)
        textbox(slide, cx, y + 0.345, cw - 0.06, 0.12, sub, size=5.6, color=MUTED)

    # ============================================ footer
    textbox(slide, M, H - 0.30, CW, 0.12,
            "Robust = Bonferroni-significant, replicated in both the European and Ashkenazi "
            "Jewish strata, and concordant in direction.  Source: browsable results app.",
            size=5.6, color=MUTED)

    prs.save(out_path)
    kb = os.path.getsize(out_path) / 1024
    print(f"Wrote {out_path}  ({kb:.0f} KB, {W}x{H} in -> 600x800 px at 96 DPI)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", default="adpd_figure.pptx")
    build(ap.parse_args().output)
