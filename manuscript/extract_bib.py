#!/usr/bin/env python
"""Extract cited BibTeX entries verbatim from a master .bib.

Reads the cited keys from a LaTeX .aux (\\citation{...}) and copies each
matching entry byte-for-byte from the master .bib into a trimmed .bib,
preserving every field (including non-standard ones like `day` that
bibexport silently drops). Entries are emitted in citation order of first
appearance so the output is deterministic.
"""
import re
import sys

aux_path, master_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

# cited keys, first-appearance order
cited = []
seen = set()
for m in re.finditer(r"\\citation\{([^}]*)\}", open(aux_path).read()):
    for key in m.group(1).split(","):
        key = key.strip()
        if key and key not in seen:
            seen.add(key)
            cited.append(key)

text = open(master_path).read()

# index every entry: @type{key, ... } with brace matching
entries = {}
for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text):
    start = m.start()
    # walk from the opening brace of the entry, matching braces
    brace_open = text.index("{", m.start())
    depth = 0
    i = brace_open
    while i < len(text):
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    key = m.group(2)
    entries[key] = text[start : i + 1]

missing = [k for k in cited if k not in entries]
if missing:
    sys.exit(f"ERROR: cited keys absent from master bib: {missing}")

# @string macros (journal abbreviations etc.), verbatim in source order,
# so any abbreviation used by a cited entry resolves.
strings = []
for m in re.finditer(r"@[Ss]tring\s*\{", text):
    brace_open = m.end() - 1
    depth = 0
    i = brace_open
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    strings.append(text[m.start() : i + 1])

with open(out_path, "w") as f:
    f.write(
        "% Cited-only bibliography for the ENSO-conditioned GMST capsule.\n"
        "% Frozen source file: this is what the manuscript builds against.\n"
        "% Extracted verbatim (every field preserved, including non-standard\n"
        "% `day` fields and @string abbreviations) from the manuscript's master\n"
        "% all.bib by manuscript/extract_bib.py. That regeneration is provenance\n"
        "% only -- it needs the external all.bib and is NOT part of `make all`.\n"
        "% Do not hand-edit.\n\n"
    )
    for s in strings:
        f.write(s.rstrip() + "\n")
    f.write("\n")
    for key in cited:
        f.write(entries[key].rstrip() + "\n\n")

print(f"wrote {len(cited)} entries to {out_path}")
