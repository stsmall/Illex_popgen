"""Export every manuscript table to editable CSV + a single multi-sheet .xlsx.

Parses the LaTeX tabular environments in main.tex, supplement.tex and software_table.tex,
strips LaTeX markup, and writes:
  tables_export/Table<N>_<label>.csv  (one per table, N as numbered in the compiled main/supplement PDFs,
                                      read from main.aux and supplement.aux; run after make all)
  tables_export/illex_tables.xlsx  (one sheet per table; caption on the first row)
Intended so the numbers/labels can be hand-edited without touching LaTeX.
"""
import os, re, csv, glob
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "tables_export")
os.makedirs(OUT, exist_ok=True)
SRC = [os.path.join(ROOT, f) for f in ("main.tex", "supplement.tex", "software_table.tex", "sampling_table.tex",
                                       "pi_blocks_table.tex", "inversion_gene_table.tex")]
NUM = {}
for aux in ("main.aux", "supplement.aux"):
    for lab, n in re.findall(r"\\newlabel\{(s?tab:[^}]+)\}\{\{(S?\d+)\}", open(os.path.join(ROOT, aux)).read()):
        NUM[lab] = n


def clean(cell):
    s = cell.strip()
    s = re.sub(r"%.*$", "", s)                                   # trailing LaTeX comment
    s = re.sub(r"\\multicolumn\{\d+\}\{[^}]*\}\{(.*?)\}", r"\1", s)
    s = re.sub(r"\\(textbf|textit|emph|texttt|textsc|mathrm|mathit|text)\{(.*?)\}", r"\2", s)
    s = re.sub(r"\\(illex|Fst|dxy|pipi)\b", lambda m: {
        "illex": "I. illecebrosus", "Fst": "FST", "dxy": "dXY", "pipi": "pi"}.get(m.group(1), ""), s)
    s = s.replace("$", "").replace("\\,", " ").replace("\\;", " ").replace("~", " ")
    s = s.replace("\\%", "%").replace("\\&", "&").replace("\\_", "_").replace("\\#", "#")
    s = s.replace("\\to", "->").replace("\\times", "x").replace("\\approx", "~")
    s = s.replace("\\alpha", "alpha").replace("\\pi", "pi").replace("\\mu", "mu")
    s = re.sub(r"\\cite[a-zA-Z]*\{[^}]*\}", "", s)              # drop citation keys entirely
    s = s.replace("\\emph", "").replace("{", "").replace("}", "")
    s = re.sub(r"\\[a-zA-Z]+", "", s)                            # any leftover macro
    s = s.replace("--", "-").replace("`", "").replace("''", '"')
    return re.sub(r"\s+", " ", s).strip()


def caption_of(block):
    m = re.search(r"\\caption\{", block)
    if not m:
        return ""
    i, depth, buf = m.end(), 1, []
    while i < len(block) and depth:
        c = block[i]
        depth += (c == "{") - (c == "}")
        if depth:
            buf.append(c)
        i += 1
    cap = clean("".join(buf))
    cap = re.sub(r"\s*\\?ref\{[^}]*\}", "", cap)
    return re.split(r"\.(?:\s|$)", cap)[0][:220]


def parse_tabular(block):
    m = re.search(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}", block, re.S) or \
        re.search(r"\\begin\{tabularx\}\{[^}]*\}\{(?:[^{}]|\{[^{}]*\})*\}(.*?)\\end\{tabularx\}", block, re.S)
    if m:
        body = m.group(1)
    else:                                   # longtable: drop caption/label and repeated head/foot
        m = re.search(r"\\begin\{longtable\}\{[^}]*\}(.*?)\\end\{longtable\}", block, re.S)
        if not m:
            return []
        body = m.group(1)
        body = body[body.index("\\endfirsthead"):] if "\\endfirsthead" in body else body
        body = re.sub(r"^.*?\\endfoot", "", body, flags=re.S) if "\\endfoot" in body else body
        body = re.sub(r"^.*?\\endhead", "", body, flags=re.S) if "\\endhead" in body else body
        hdr = re.search(r"\\toprule\s*(.*?)\\\\\s*\\midrule", m.group(1), re.S)
        body = (hdr.group(1) + r" \\ " if hdr else "") + body
    body = re.sub(r"\\(toprule|midrule|bottomrule|hline|addlinespace)(\[[^]]*\])?", "", body)
    body = re.sub(r"\\cmidrule(\([^)]*\))?\{[^}]*\}", "", body)
    rows = []
    for raw in body.split(r"\\"):
        raw = raw.strip()
        if not raw:
            continue
        cells = [clean(c) for c in re.split(r"(?<!\\)&", raw)]
        if any(cells):
            rows.append(cells)
    return rows


def label_of(block, i):
    m = re.search(r"\\label\{(s?tab:[^}]*)\}", block)
    if not m:
        return f"table{i}"
    return f"Table{NUM[m.group(1)]}_{m.group(1).split(':')[1]}"


tables = []
for path in SRC:
    if not os.path.exists(path):
        continue
    txt = open(path).read()
    for block in re.findall(r"\\begin\{(?:long)?table\*?\}.*?\\end\{(?:long)?table\*?\}", txt, re.S):
        rows = parse_tabular(block)
        if not rows:
            continue
        tables.append((label_of(block, len(tables) + 1), caption_of(block), rows))
    # software_table.tex may hold a bare tabular with no table env
    if "software_table" in path and not re.search(r"\\begin\{table", txt):
        rows = parse_tabular(txt)
        if rows:
            tables.append((f"Table{NUM['tab:software']}_software", "Software and tools used", rows))

def order(t):
    n = t[0][5:].split("_")[0]
    return (n.startswith("S"), int(n.lstrip("S")) if n.lstrip("S").isdigit() else 99)
tables.sort(key=order)

# ---- write CSVs ----
for label, cap, rows in tables:
    with open(os.path.join(OUT, f"{label}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        if cap:
            w.writerow([cap])
        w.writerows(rows)

# ---- write one xlsx workbook ----
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
wb = Workbook()
wb.remove(wb.active)
for label, cap, rows in tables:
    ws = wb.create_sheet(title=label[:31])
    if cap:
        ws.append([cap])
        ws["A1"].font = Font(italic=True, size=9)
        ws["A1"].alignment = Alignment(wrap_text=False)
    hdr_row = ws.max_row + 1
    for r in rows:
        ws.append(r)
    for c in ws[hdr_row]:            # bold the header row (first tabular row)
        c.font = Font(bold=True)
    for col in ws.columns:
        width = max((len(str(c.value)) for c in col if c.value), default=10)
        ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)
wb.save(os.path.join(OUT, "illex_tables.xlsx"))

print(f"wrote {len(tables)} tables to {OUT}")
for label, cap, rows in tables:
    print(f"  {label:14s} {len(rows):2d} rows x {max(len(r) for r in rows)} cols  | {cap[:60]}")
