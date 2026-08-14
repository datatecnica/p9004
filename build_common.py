"""Shared helpers for the Project 9004 rebuild.

The important thing in here is `make_key()`. Every merge in this pipeline joins on
`key = PATNO + '_' + VISIT`, and the failure mode we are guarding against is a silent
one: a stray space or an Excel float turns `####_BL` into `#### _BL` or `####.0_BL`,
the left join misses, and the block arrives 100% missing without anything erroring.

So `make_key` normalises, and records every normalisation it had to perform. Defects
that are unambiguous typos (whitespace, `.0` suffix) are repaired and reported.
Defects that would require guessing (a non-numeric PATNO, an empty visit) are raised.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime

import pandas as pd

# Characters that survive a copy-paste through Excel and are invisible in a terminal.
_INVISIBLE = {
    " ": "NBSP",           # non-breaking space
    "​": "ZWSP",           # zero-width space
    "‌": "ZWNJ",
    "‍": "ZWJ",
    "﻿": "BOM",
    " ": "FIGSP",
    " ": "NNBSP",
}

_PATNO_OK = re.compile(r"^\d+$")
_VISIT_OK = re.compile(r"^[A-Z0-9]+$")
_FLOATY = re.compile(r"^(\d+)\.0+$")

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def log(msg: str = "") -> None:
    print(msg, flush=True)


def _scrub(s: pd.Series, findings: list, source: str, col: str) -> pd.Series:
    """Normalise one key component, appending a row to `findings` per defect class."""
    out = s.astype("string")

    # Missing before we do anything else.
    blank = out.isna() | out.str.strip().isin(["", "nan", "NaN", "None", "NA", "<NA>"])
    if blank.any():
        findings.append(
            dict(source=source, column=col, defect="empty_or_null",
                 n=int(blank.sum()), action="FATAL",
                 example=repr(out[blank].iloc[0]) if out[blank].notna().any() else "<NA>")
        )

    # Invisible characters.
    for ch, name in _INVISIBLE.items():
        hit = out.str.contains(ch, regex=False, na=False)
        if hit.any():
            findings.append(
                dict(source=source, column=col, defect=f"invisible_char_{name}",
                     n=int(hit.sum()), action="repaired",
                     example=repr(out[hit].iloc[0]))
            )
            out = out.str.replace(ch, "", regex=False)

    # Leading / trailing whitespace.
    ws = out.notna() & (out != out.str.strip())
    if ws.any():
        findings.append(
            dict(source=source, column=col, defect="surrounding_whitespace",
                 n=int(ws.sum()), action="repaired", example=repr(out[ws].iloc[0]))
        )
    out = out.str.strip()
    return out


def make_key(
    df: pd.DataFrame,
    source: str,
    patno_col: str = "PATNO",
    visit_col: str = "EVENT_ID",
) -> tuple[pd.Series, pd.DataFrame]:
    """Build `key` = PATNO_VISIT, returning (key, findings).

    Repairs whitespace/invisible-char/float-formatting defects and reports them.
    Raises on defects that cannot be repaired without guessing.
    """
    findings: list = []

    patno = _scrub(df[patno_col], findings, source, patno_col)
    visit = _scrub(df[visit_col], findings, source, visit_col)

    # Excel/pandas float formatting: 3000.0 -> 3000
    floaty = patno.str.match(_FLOATY, na=False)
    if floaty.any():
        findings.append(
            dict(source=source, column=patno_col, defect="float_formatted_id",
                 n=int(floaty.sum()), action="repaired",
                 example=repr(patno[floaty].iloc[0]))
        )
        patno = patno.str.replace(_FLOATY, r"\1", regex=True)

    # Visit-code case drift (v04 vs V04). Uppercase is the PPMI convention.
    lower = visit.notna() & (visit != visit.str.upper())
    if lower.any():
        findings.append(
            dict(source=source, column=visit_col, defect="lowercase_visit_code",
                 n=int(lower.sum()), action="repaired", example=repr(visit[lower].iloc[0]))
        )
    visit = visit.str.upper()

    # Anything still malformed is not a typo we can fix.
    bad_p = patno.notna() & ~patno.str.match(_PATNO_OK, na=False)
    if bad_p.any():
        findings.append(
            dict(source=source, column=patno_col, defect="non_numeric_id",
                 n=int(bad_p.sum()), action="FATAL", example=repr(patno[bad_p].iloc[0]))
        )
    bad_v = visit.notna() & ~visit.str.match(_VISIT_OK, na=False)
    if bad_v.any():
        findings.append(
            dict(source=source, column=visit_col, defect="malformed_visit_code",
                 n=int(bad_v.sum()), action="FATAL", example=repr(visit[bad_v].iloc[0]))
        )

    key = (patno + "_" + visit).rename("key")

    dup = key.duplicated(keep=False) & key.notna()
    if dup.any():
        findings.append(
            dict(source=source, column="key", defect="duplicate_key",
                 n=int(dup.sum()), action="FATAL", example=repr(key[dup].iloc[0]))
        )

    fdf = pd.DataFrame(findings, columns=["source", "column", "defect", "n", "action", "example"])
    return key, fdf


def make_patno(
    df: pd.DataFrame,
    source: str,
    patno_col: str = "PATNO",
) -> tuple[pd.Series, pd.DataFrame]:
    """PATNO-only variant of `make_key`, for participant-level blocks (genetics).

    Same normalisation rules; duplicates are fatal here too, since a participant-level
    block with a repeated PATNO would fan out the scaffold on join.
    """
    findings: list = []
    patno = _scrub(df[patno_col], findings, source, patno_col)

    floaty = patno.str.match(_FLOATY, na=False)
    if floaty.any():
        findings.append(
            dict(source=source, column=patno_col, defect="float_formatted_id",
                 n=int(floaty.sum()), action="repaired",
                 example=repr(patno[floaty].iloc[0]))
        )
        patno = patno.str.replace(_FLOATY, r"\1", regex=True)

    bad = patno.notna() & ~patno.str.match(_PATNO_OK, na=False)
    if bad.any():
        findings.append(
            dict(source=source, column=patno_col, defect="non_numeric_id",
                 n=int(bad.sum()), action="FATAL", example=repr(patno[bad].iloc[0]))
        )

    dup = patno.duplicated(keep=False) & patno.notna()
    if dup.any():
        findings.append(
            dict(source=source, column=patno_col, defect="duplicate_patno",
                 n=int(dup.sum()), action="FATAL", example=repr(patno[dup].iloc[0]))
        )

    fdf = pd.DataFrame(findings, columns=["source", "column", "defect", "n", "action", "example"])
    return patno.rename("PATNO"), fdf


def assert_no_fatal(findings: pd.DataFrame) -> None:
    """Stop the build on any defect that cannot be repaired by rule."""
    if findings.empty:
        return
    fatal = findings[findings["action"] == "FATAL"]
    if not fatal.empty:
        log("\nFATAL key defects — refusing to merge:")
        log(fatal.to_string(index=False))
        sys.exit(1)


def report_findings(findings: pd.DataFrame, label: str) -> None:
    if findings.empty:
        log(f"  key hygiene [{label}]: clean")
    else:
        log(f"  key hygiene [{label}]: {len(findings)} defect class(es)")
        for _, r in findings.iterrows():
            log(f"    {r['action']:9s} {r['defect']:26s} n={r['n']:<7d} e.g. {r['example']}")
