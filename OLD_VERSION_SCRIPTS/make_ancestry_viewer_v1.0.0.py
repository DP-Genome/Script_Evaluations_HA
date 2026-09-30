#!/usr/bin/env python3
"""
make_ancestry_viewer.py  (version 1.0.0)

Builds a single, self-contained HTML page that shows the PCA and ADMIXTURE results of one
sample against the HGDP reference panel. It needs only the Python 3 standard library and
reads the files written by pileupCaller_glimpse.sh / pca_admixture.sh (v1.3.0 or later):

    <run_dir>/smart_pca/merged_ref<suffix>_labeled.evec   PCA coordinates, one row per individual
    <run_dir>/smart_pca/merged_ref<suffix>.eval           eigenvalues (percent variance per PC)
    <run_dir>/confEllipses                                95% ellipse of the sample (optional)
    <run_dir>/admixture<suffix>/Qfiles_withpop/*.Q.withPop  annotated Q files, one per K (optional)
    <run_dir>/admixture<suffix>/cv_comparison.txt          CV errors (optional)
    <run_dir>/pipeline_*.log                               SNP counts and coverage (optional)

<suffix> is "_tv" for runs made with --transversions true, otherwise empty. It is detected
automatically and can be forced with --suffix.

Examples
    python3 make_ancestry_viewer.py /path/to/expo_250926
    python3 make_ancestry_viewer.py . --suffix "" --case "Casework 2026" --out viewer.html

The page has no external data: open the .html file in any browser. Fonts come from Google
Fonts when online and fall back to system fonts otherwise.
"""
import argparse
import glob
import html
import json
import os
import re
import sys

VERSION = "1.0.0"
# fixed continent order and display names; the page assigns colours and marker shapes by position
CONT = ["Africa", "Europe", "Middle_East", "Central_South_Asia", "East_Asia", "Oceania", "America"]
LABEL = {"Africa": "Africa", "Europe": "Europe", "Middle_East": "Middle East",
         "Central_South_Asia": "Central & South Asia", "East_Asia": "East Asia",
         "Oceania": "Oceania", "America": "America"}
CI = {c: i for i, c in enumerate(CONT)}


def die(msg):
    sys.exit("make_ancestry_viewer: " + msg)


def to_int(s):
    return int(s.replace(",", ""))


def find_suffix(run, forced):
    if forced is not None:
        return forced
    for suffix in ("_tv", ""):
        if os.path.exists(f"{run}/smart_pca/merged_ref{suffix}_labeled.evec"):
            return suffix
    die(f"no smart_pca/merged_ref*_labeled.evec found in {run}. Run the pipeline with --plot pca or both first.")


def read_pca(run, suffix, want):
    evec = f"{run}/smart_pca/merged_ref{suffix}_labeled.evec"
    evalf = f"{run}/smart_pca/merged_ref{suffix}.eval"
    if not os.path.exists(evec):
        die(f"missing {evec}")
    if not os.path.exists(evalf):
        die(f"missing {evalf}")
    ev = [float(x) for x in open(evalf).read().split()]
    tot = sum(ev)
    subs, pts, samples, npc = [], [], [], None
    for line in open(evec):
        if line.lstrip().startswith("#") or not line.strip():
            continue
        parts = line.split()
        # ID PC1..PCn population subpopulation
        nums = []
        for tok in parts[1:]:
            try:
                nums.append(float(tok))
            except ValueError:
                break
        n = len(nums)
        if npc is None:
            npc = n
        cont = parts[1 + n] if len(parts) > 1 + n else "?"
        sub = parts[2 + n] if len(parts) > 2 + n else parts[0]
        pid = parts[0]
        if cont in CI:
            if sub not in subs:
                subs.append(sub)
            pts.append([CI[cont], subs.index(sub), pid.split(":")[-1]] + [round(x, 4) for x in nums])
        else:
            samples.append({"id": pid.split(":")[-1], "full": pid, "pcs": nums})
    if not samples:
        die("no non-reference individual found in the evec file")
    ids = [s["id"] for s in samples]
    if want:
        pick = [s for s in samples if s["id"] == want]
        if not pick:
            die(f"--sample {want} not found. Samples in the evec file: {', '.join(ids)}")
        sample = pick[0]
    elif len(samples) == 1:
        sample = samples[0]
    else:
        die(f"several non-reference samples found ({', '.join(ids)}). Choose one with --sample.")
    pct = [round(100 * e / tot, 2) for e in ev[:npc]]
    return {"pct": pct, "subs": subs, "pts": pts, "sample": sample}


def read_ellipse(run, sample_full):
    path = f"{run}/confEllipses"
    if not os.path.exists(path):
        return None
    current, found = None, None
    for line in open(path):
        if line.startswith("sample:"):
            current = line.split(":", 1)[1].strip()
        elif line.startswith("ellcoords:") and current == sample_full:
            v = line.split(":", 1)[1].split("::")[0].split()
            conf = float(line.rsplit("::", 1)[1])
            found = {"cx": float(v[0]), "cy": float(v[1]), "a": float(v[2]), "b": float(v[3]),
                     "ang": float(v[4]), "conf": conf}
    return found


def read_admixture(run, suffix, sample):
    adir = f"{run}/admixture{suffix}"
    files = sorted(glob.glob(f"{adir}/Qfiles_withpop/*.Q.withPop"))
    if not files:
        return None
    cv = {}
    cvfile = f"{adir}/cv_comparison.txt"
    if os.path.exists(cvfile):
        for line in open(cvfile):
            m = re.search(r"K=(\d+)\):\s*([\d.]+)", line)
            if m:
                cv[int(m.group(1))] = float(m.group(2))
    upper = {c.upper(): c for c in CONT}
    pops, byk, counts = {}, {}, {}
    for f in files:
        m = re.search(r"\.(\d+)\.Q\.withPop$", f)
        if not m:
            continue
        k = int(m.group(1))
        acc, smp = {}, None
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            q = [float(x) for x in " ".join(p[4:]).split()]
            if p[2] in upper:
                pops[p[1]] = upper[p[2]]
                a = acc.setdefault(p[1], [0, [0.0] * k])
                a[0] += 1
                for i, x in enumerate(q):
                    a[1][i] += x
            elif p[0].split(":")[-1] == sample:
                smp = q
        if smp is None:
            die(f"sample {sample} not found in {f}")
        byk[k] = {"sample": [round(x, 4) for x in smp],
                  "avg": {n: [round(x / a[0], 4) for x in a[1]] for n, a in acc.items()}}
        counts = {n: a[0] for n, a in acc.items()}
    names = sorted(pops, key=lambda n: (CI[pops[n]], n.lower()))
    for k in byk:
        byk[k]["avg"] = [byk[k]["avg"][n] for n in names]
    for k in byk:
        cv.setdefault(k, None)
    if any(v is None for v in cv.values()):
        die(f"{cvfile} is missing or has no CV error for every K")
    return {"cv": cv, "pops": [{"name": n, "cont": CI[pops[n]], "n": counts[n]} for n in names], "byK": byk}


def read_meta(run, suffix, sample, log, n_individuals):
    meta = {"individuals": n_individuals, "snps": None, "covered": None, "coverage": None}
    if log is None:
        logs = [f for f in glob.glob(f"{run}/pipeline_*.log") if "aborted" not in f]
        log = max(logs, key=os.path.getmtime) if logs else None
    text = open(log, errors="replace").read() if log and os.path.exists(log) else ""
    m = re.findall(r"final: ([\d,]+) SNPs[^\n]*?priority retained ([\d,]+)/", text)
    if m:
        meta["snps"], meta["covered"] = to_int(m[-1][0]), to_int(m[-1][1])
    else:
        bim = f"{run}/merged_plink_ldprio/merged_ref{suffix}_ldprioPruned.bim"
        if os.path.exists(bim):
            meta["snps"] = sum(1 for _ in open(bim))
    m = re.search(r"^" + re.escape(sample) + r"\t\d+\t\d+\t([0-9.eE+-]+)", text, re.M)
    if m:
        meta["coverage"] = round(float(m.group(1)), 4)
    if suffix == "_tv":
        m = re.search(r"Transversions kept: (\d+)", text)
        if m:
            meta["transversionsKept"] = int(m.group(1))
    return meta


def main():
    ap = argparse.ArgumentParser(description="Build a self-contained PCA + ADMIXTURE viewer page.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Examples")[0])
    ap.add_argument("run_dir", nargs="?", default=".", help="pipeline run directory (default: current directory)")
    ap.add_argument("--sample", help="sample ID to highlight (default: the only non-reference sample)")
    ap.add_argument("--suffix", help='"_tv" for transversion runs, "" otherwise (default: detect)')
    ap.add_argument("--case", help="short case label for the page header (default: run directory name)")
    ap.add_argument("--log", help="pipeline log for SNP counts and coverage (default: newest pipeline_*.log)")
    ap.add_argument("--out", help="output HTML file (default: <run_dir>/<sample>_ancestry_viewer.html)")
    ap.add_argument("--version", action="version", version=VERSION)
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    if not os.path.isdir(run):
        die(f"{run} is not a directory")
    suffix = find_suffix(run, args.suffix)
    pca = read_pca(run, suffix, args.sample)
    sample = pca["sample"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample["id"]):
        die(f"sample ID '{sample['id']}' has characters other than letters, digits, . _ -")
    pca["ellipse"] = read_ellipse(run, sample["full"])
    admix = read_admixture(run, suffix, sample["id"])
    meta = read_meta(run, suffix, sample["id"], args.log, len(pca["pts"]) + 1)
    case = args.case or os.path.basename(run).replace("_", " ")
    data = {"meta": meta,
            "continents": [{"key": c, "label": LABEL[c], "n": sum(1 for p in pca["pts"] if p[0] == CI[c])} for c in CONT],
            "pca": pca, "admix": admix}
    page = (TEMPLATE.replace("__SAMPLE__", sample["id"]).replace("__CASE__", html.escape(case))
            .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    out = args.out or f"{run}/{sample['id']}_ancestry_viewer.html"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(page)
    print(f"sample {sample['id']} · {len(pca['pts'])} reference individuals · {len(pca['subs'])} populations · "
          f"{len(pca['sample']['pcs'])} PCs · ellipse: {'yes' if pca['ellipse'] else 'no'} · "
          f"ADMIXTURE K={','.join(map(str, sorted(admix['byK']))) if admix else 'none'}")
    print(f"wrote {out} ({os.path.getsize(out) // 1024} KB)")


TEMPLATE = r'''<title>__SAMPLE__ Ancestry</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Instrument+Sans:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
/* Layout: an instrument panel. Slim header with the mode switch, a controls rail on the left, the plot or bars on the right.
   Palette: cool paper greys leaning navy, ink for text and for the case sample (so it never competes with a population hue),
   seven validated categorical hues for the HGDP continents. Type: Instrument Sans for the interface, JetBrains Mono for numbers. */
:root {
  --paper: #e9edf3;
  --surface: #fbfcfd;
  --raised: #f2f5f9;
  --ink: #0e1a2c;
  --ink-2: #46546a;
  --ink-3: #75839a;
  --line: #d3dae5;
  --grid: #e3e8ef;
  --focus: #0e1a2c;
  --c0: #2a78d6; /* Africa */
  --c1: #e87ba4; /* Europe */
  --c2: #008300; /* Middle East */
  --c3: #eda100; /* Central & South Asia */
  --c4: #4a3aa7; /* East Asia */
  --c5: #eb6834; /* Oceania */
  --c6: #1baf7a; /* America */
  --font-ui: "Instrument Sans", "Segoe UI", system-ui, -apple-system, sans-serif;
  --font-data: "JetBrains Mono", ui-monospace, "SF Mono", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --paper: #0c1118;
    --surface: #141b25;
    --raised: #1a2330;
    --ink: #e8edf5;
    --ink-2: #aab6c8;
    --ink-3: #7d8aa0;
    --line: #2a3546;
    --grid: #1f2937;
    --focus: #e8edf5;
    --c0: #3987e5; --c1: #d55181; --c2: #008300; --c3: #c98500; --c4: #9085e9; --c5: #d95926; --c6: #199e70;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --paper: #0c1118;
  --surface: #141b25;
  --raised: #1a2330;
  --ink: #e8edf5;
  --ink-2: #aab6c8;
  --ink-3: #7d8aa0;
  --line: #2a3546;
  --grid: #1f2937;
  --focus: #e8edf5;
  --c0: #3987e5; --c1: #d55181; --c2: #008300; --c3: #c98500; --c4: #9085e9; --c5: #d95926; --c6: #199e70;
}

* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--paper);
  color: var(--ink);
  font: 400 14px/1.45 var(--font-ui);
  padding-inline: 20px;
}
.app { max-width: 1280px; margin-inline: auto; padding-block: 22px 40px; display: flex; flex-direction: column; gap: 18px; }
button, input, select { font: inherit; color: inherit; }
button { cursor: pointer; }
:focus-visible { outline: 2px solid var(--focus); outline-offset: 2px; }
.mono { font-family: var(--font-data); font-variant-numeric: tabular-nums; }

/* header */
.head { display: flex; flex-wrap: wrap; align-items: flex-end; justify-content: space-between; gap: 16px 32px; }
.eyebrow { font: 500 11px/1 var(--font-data); letter-spacing: .09em; text-transform: uppercase; color: var(--ink-3); margin: 0 0 8px; }
h1 { margin: 0; font-size: 30px; line-height: 1.05; font-weight: 700; letter-spacing: -.02em; text-wrap: balance; }
.facts { display: flex; flex-wrap: wrap; gap: 6px 22px; margin: 12px 0 0; padding: 0; list-style: none; color: var(--ink-2); font-size: 13px; }
.facts b { font-family: var(--font-data); font-weight: 500; color: var(--ink); }
.modes { display: inline-flex; border: 1px solid var(--ink); border-radius: 4px; overflow: hidden; }
.modes button { border: 0; background: transparent; padding: 9px 22px; font-weight: 600; letter-spacing: .04em; font-size: 13px; }
.modes button[aria-pressed="true"] { background: var(--ink); color: var(--surface); }

/* panels */
.panel { background: var(--surface); border: 1px solid var(--line); border-radius: 6px; }
.work { display: grid; grid-template-columns: 288px minmax(0, 1fr); gap: 18px; align-items: start; }
.rail { display: flex; flex-direction: column; gap: 14px; padding: 16px; }
.rail h2 { margin: 0 0 8px; font: 500 11px/1 var(--font-data); letter-spacing: .09em; text-transform: uppercase; color: var(--ink-3); }
.rail section + section { border-top: 1px solid var(--line); padding-top: 14px; }
.stage { padding: 14px 14px 10px; min-width: 0; position: relative; }

/* segmented buttons */
.pcgrid { display: grid; grid-template-columns: repeat(5, 1fr); gap: 4px; }
.chip { border: 1px solid var(--line); background: var(--surface); border-radius: 4px; padding: 6px 0; font: 500 12px/1 var(--font-data); text-align: center; }
.chip:hover { border-color: var(--ink-3); }
.chip[aria-pressed="true"] { background: var(--ink); border-color: var(--ink); color: var(--surface); }
.row-btns { display: flex; flex-wrap: wrap; gap: 6px; }
.btn { border: 1px solid var(--line); background: var(--surface); border-radius: 4px; padding: 7px 11px; font-weight: 500; font-size: 13px; }
.btn:hover { border-color: var(--ink-3); background: var(--raised); }
.btn:disabled { color: var(--ink-3); cursor: not-allowed; background: var(--surface); }
.btn.icon { width: 34px; padding-inline: 0; font-family: var(--font-data); font-size: 16px; line-height: 1; }
.toggle { display: flex; align-items: flex-start; gap: 10px; }
.toggle input { margin: 3px 0 0; width: 16px; height: 16px; accent-color: var(--ink); flex: none; }
.toggle label { font-weight: 500; }
.hint { color: var(--ink-3); font-size: 12px; margin: 4px 0 0; }

/* continent legend */
.legend { display: flex; flex-direction: column; gap: 2px; }
.legend button { display: flex; align-items: center; gap: 10px; width: 100%; border: 0; background: transparent; padding: 5px 6px; border-radius: 4px; text-align: left; }
.legend button:hover { background: var(--raised); }
.legend button[aria-pressed="false"] { opacity: .38; }
.legend .n { margin-left: auto; color: var(--ink-3); font: 400 12px var(--font-data); }
.legend svg { flex: none; }

/* PCA plot */
.plotwrap { position: relative; }
#pca { display: block; width: 100%; touch-action: none; user-select: none; -webkit-user-select: none; cursor: crosshair; }
.axis-label { font: 500 12px var(--font-ui); fill: var(--ink-2); }
.tick { font: 400 11px var(--font-data); fill: var(--ink-3); }
.tip { position: absolute; pointer-events: none; background: var(--ink); color: var(--surface); padding: 8px 10px; border-radius: 4px; font-size: 12px; line-height: 1.4; max-width: 240px; z-index: 5; }
.tip b { display: block; font-size: 13px; font-weight: 600; }
.tip .m { font-family: var(--font-data); font-variant-numeric: tabular-nums; }
.tip .k { display: flex; align-items: center; gap: 7px; }
.tip .k i { width: 12px; height: 8px; border-radius: 1px; flex: none; outline: 1px solid color-mix(in srgb, var(--surface) 40%, transparent); }
.toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; margin-top: 10px; }
.toolbar .hint { margin: 0 0 0 auto; }
.readout { margin: 0; display: grid; grid-template-columns: auto 1fr; gap: 3px 12px; font-size: 13px; }
.readout dt { color: var(--ink-2); }
.readout dd { margin: 0; font-family: var(--font-data); font-variant-numeric: tabular-nums; }
.near { margin: 8px 0 0; padding: 0; list-style: none; display: flex; flex-direction: column; gap: 2px; font-size: 13px; }
.near li { display: flex; gap: 8px; align-items: center; }
.near .d { margin-left: auto; color: var(--ink-3); font: 400 12px var(--font-data); }

/* ADMIXTURE */
.ks { display: flex; flex-wrap: wrap; gap: 6px; }
.kbtn { display: flex; flex-direction: column; align-items: center; gap: 3px; min-width: 66px; padding: 8px 10px 7px; border: 1px solid var(--line); background: var(--surface); border-radius: 4px; }
.kbtn:hover { border-color: var(--ink-3); }
.kbtn b { font: 600 15px/1 var(--font-data); }
.kbtn small { font: 400 10.5px/1 var(--font-data); color: var(--ink-3); }
.kbtn[aria-pressed="true"] { background: var(--ink); border-color: var(--ink); color: var(--surface); }
.kbtn[aria-pressed="true"] small { color: color-mix(in srgb, var(--surface) 70%, transparent); }
.kbtn .low { font: 600 9.5px/1 var(--font-data); letter-spacing: .06em; text-transform: uppercase; color: var(--c2); }
.kbtn[aria-pressed="true"] .low { color: var(--surface); }
.summary { display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 18px; margin: 14px 0 4px; }
.summary h3 { margin: 0; font-size: 15px; font-weight: 600; }
.summary .part { display: inline-flex; align-items: center; gap: 6px; font-size: 13px; }
.summary .part i { width: 10px; height: 10px; border-radius: 2px; }
.comp-legend { display: flex; flex-wrap: wrap; gap: 4px 16px; margin: 8px 0 12px; font-size: 12px; color: var(--ink-2); }
.comp-legend span { display: inline-flex; align-items: center; gap: 6px; }
.comp-legend i { width: 18px; height: 10px; border-radius: 2px; outline: 1px solid var(--line); outline-offset: -1px; }
.view-tabs { display: inline-flex; border: 1px solid var(--line); border-radius: 4px; overflow: hidden; margin-left: auto; }
.view-tabs button { border: 0; background: var(--surface); padding: 5px 12px; font-size: 12.5px; font-weight: 500; }
.view-tabs button[aria-pressed="true"] { background: var(--ink); color: var(--surface); }
.stage-head { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 12px; }

.bars { display: flex; flex-direction: column; }
.grp { font: 500 10.5px/1 var(--font-data); letter-spacing: .09em; text-transform: uppercase; color: var(--ink-3); margin: 12px 0 5px 158px; }
.brow { display: flex; align-items: center; gap: 8px; padding-block: 1px; }
.brow .lab { width: 150px; flex: none; text-align: right; font-size: 12.5px; color: var(--ink-2); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.brow .lab small { color: var(--ink-3); font: 400 10.5px var(--font-data); margin-left: 4px; }
.bar { flex: 1; display: flex; height: 17px; min-width: 0; background: var(--raised); border-radius: 2px; overflow: hidden; }
.bar i { display: block; height: 100%; flex: none; box-shadow: inset -2px 0 0 var(--surface); }
.bar i:hover { filter: brightness(1.12) saturate(1.1); }
.brow.me { padding-block: 4px 8px; margin-bottom: 4px; border-bottom: 1px solid var(--line); }
.brow.me .lab { font-weight: 700; color: var(--ink); font-size: 14px; }
.brow.me .bar { height: 30px; outline: 2px solid var(--ink); outline-offset: 1px; }
.empty { padding: 40px 0; text-align: center; color: var(--ink-3); }

.picker input[type="search"] { width: 100%; padding: 7px 9px; border: 1px solid var(--line); border-radius: 4px; background: var(--surface); margin-bottom: 8px; }
.plist { max-height: calc(100vh - 230px); min-height: 220px; overflow: auto; padding-right: 4px; }
.picker { position: sticky; top: 12px; }
.pgroup { margin-bottom: 8px; }
.pgroup > .gh { display: flex; align-items: center; gap: 8px; font-weight: 600; font-size: 13px; padding: 4px 0; }
.pgroup label.p { display: flex; align-items: center; gap: 8px; padding: 2px 0 2px 22px; font-size: 13px; }
.pgroup label.p .n { margin-left: auto; color: var(--ink-3); font: 400 11.5px var(--font-data); }
.pgroup input { accent-color: var(--ink); width: 15px; height: 15px; margin: 0; flex: none; }
.sw { width: 10px; height: 10px; border-radius: 2px; flex: none; }

.tablewrap { overflow-x: auto; border: 1px solid var(--line); border-radius: 4px; }
table { border-collapse: collapse; width: 100%; font-size: 12.5px; }
th, td { padding: 5px 9px; text-align: right; white-space: nowrap; border-bottom: 1px solid var(--line); }
th:first-child, td:first-child { text-align: left; position: sticky; left: 0; background: var(--surface); }
thead th { font-weight: 500; color: var(--ink-2); vertical-align: bottom; }
thead th i { display: block; height: 6px; margin-bottom: 4px; border-radius: 1px; }
td { font-family: var(--font-data); font-variant-numeric: tabular-nums; }
tr.me td, tr.me th { font-weight: 700; background: var(--raised); }
td.dim { color: var(--ink-3); }
.note { color: var(--ink-2); font-size: 12.5px; margin: 12px 0 0; max-width: 78ch; }

@media (max-width: 900px) {
  .work { grid-template-columns: minmax(0, 1fr); }
  .plist { max-height: 300px; }
  .picker { position: static; }
  .brow .lab { width: 104px; font-size: 11.5px; }
  .grp { margin-left: 112px; }
  .toolbar .hint { margin-left: 0; }
  h1 { font-size: 25px; }
}
</style>

<div class="app">
  <header class="head">
    <div>
      <p class="eyebrow">__CASE__ · HGDP reference panel</p>
      <h1>__SAMPLE__ ancestry</h1>
      <ul class="facts" id="facts"></ul>
    </div>
    <div class="modes" role="group" aria-label="Analysis">
      <button id="mode-pca" type="button" aria-pressed="true">PCA</button>
      <button id="mode-adm" type="button" aria-pressed="false">ADMIXTURE</button>
    </div>
  </header>

  <!-- PCA -->
  <div class="work" id="view-pca">
    <aside class="panel rail">
      <section>
        <h2>X axis</h2>
        <div class="pcgrid" id="xpcs"></div>
      </section>
      <section>
        <h2>Y axis</h2>
        <div class="pcgrid" id="ypcs"></div>
        <div class="row-btns" style="margin-top:8px"><button class="btn" id="swap" type="button">Swap axes</button></div>
      </section>
      <section>
        <div class="toggle">
          <input type="checkbox" id="ell" checked>
          <div><label for="ell">95% confidence ellipse</label><p class="hint" id="ellhint"></p></div>
        </div>
      </section>
      <section>
        <h2>Reference populations</h2>
        <div class="legend" id="legend"></div>
        <div class="row-btns" style="margin-top:8px"><button class="btn" id="showall" type="button">Show all</button></div>
      </section>
      <section>
        <h2>__SAMPLE__ in this view</h2>
        <dl class="readout" id="readout"></dl>
        <ul class="near" id="near"></ul>
      </section>
    </aside>
    <main class="panel stage">
      <div class="plotwrap" id="plotwrap"><svg id="pca" role="img" aria-label="Principal component plot of __SAMPLE__ against HGDP reference individuals"></svg><div class="tip" id="tip" hidden></div></div>
      <div class="toolbar">
        <button class="btn icon" id="zin" type="button" aria-label="Zoom in">+</button>
        <button class="btn icon" id="zout" type="button" aria-label="Zoom out">−</button>
        <button class="btn" id="zsample" type="button">Zoom to __SAMPLE__</button>
        <button class="btn" id="zreset" type="button">Reset view</button>
        <p class="hint">Drag to zoom into a region · click to zoom in 2× · Shift+drag to pan</p>
      </div>
    </main>
  </div>

  <!-- ADMIXTURE -->
  <div class="work" id="view-adm" hidden>
    <aside class="panel rail picker">
      <section>
        <h2>Populations shown</h2>
        <input type="search" id="psearch" placeholder="Search populations" aria-label="Search populations">
        <div class="row-btns" style="margin-bottom:10px">
          <button class="btn" id="pall" type="button">All</button>
          <button class="btn" id="pnone" type="button">None</button>
        </div>
        <div class="plist" id="plist"></div>
      </section>
    </aside>
    <main class="panel stage">
      <div class="stage-head">
        <div class="ks" id="ks" role="group" aria-label="Number of ancestral components (K)"></div>
        <div class="view-tabs" role="group" aria-label="View"><button id="v-bars" type="button" aria-pressed="true">Bars</button><button id="v-table" type="button" aria-pressed="false">Table</button></div>
      </div>
      <div class="summary" id="summary"></div>
      <div class="comp-legend" id="comps"></div>
      <div id="admbody"></div>
      <p class="note" id="admnote"></p>
      <div class="tip" id="tip2" hidden></div>
    </main>
  </div>
</div>

<script>
const DATA = __DATA__;
</script>
<script>
(function () {
  "use strict";
  var D = DATA, CONT = D.continents, NC = CONT.length;
  var $ = function (id) { return document.getElementById(id); };
  function el(tag, cls, text) { var e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; }
  function pct(v, d) { return (v * 100).toFixed(d == null ? 1 : d) + "%"; }

  /* ---------- header facts ---------- */
  var m = D.meta;
  var facts = $("facts");
  function fact(parts) { var li = el("li"); parts.forEach(function (p) { if (p.b) li.appendChild(el("b", null, p.b)); else li.appendChild(document.createTextNode(p)); }); facts.appendChild(li); }
  fact([{ b: m.individuals.toLocaleString() }, " individuals"]);
  fact([{ b: m.snps.toLocaleString() }, " SNPs, transversions only"]);
  (function () {
    var parts = ["__SAMPLE__ "];
    if (m.coverage != null) parts.push({ b: m.coverage.toFixed(2) + "×" }, " coverage");
    if (m.covered != null) parts.push(m.coverage != null ? ", " : "", { b: m.covered.toLocaleString() }, " SNPs called");
    if (parts.length > 1) fact(parts);
  })();

  if (!D.admix) $("mode-adm").hidden = true;

  /* ---------- mode switch ---------- */
  var mode = "pca";
  function setMode(x) {
    mode = x;
    $("mode-pca").setAttribute("aria-pressed", x === "pca");
    $("mode-adm").setAttribute("aria-pressed", x === "adm");
    $("view-pca").hidden = x !== "pca";
    $("view-adm").hidden = x !== "adm";
    if (x === "pca") drawPCA(); else drawAdm();
  }
  $("mode-pca").onclick = function () { setMode("pca"); };
  $("mode-adm").onclick = function () { setMode("adm"); };

  /* ================= PCA ================= */
  var P = D.pca, pts = P.pts, S = P.sample, NPC = S.pcs.length;
  var SHAPES = ["circle", "square", "diamond", "tri", "trid", "hex", "plus"];
  var px = 0, py = 1, showEll = true, hidden = {}, view = null, base = null;

  function shapePath(kind, x, y, r) {
    var f = function (n) { return n.toFixed(1); };
    switch (kind) {
      case "circle": return "M" + f(x - r) + "," + f(y) + "a" + r + "," + r + " 0 1,0 " + f(2 * r) + ",0a" + r + "," + r + " 0 1,0 " + f(-2 * r) + ",0z";
      case "square": return "M" + f(x - r * .9) + "," + f(y - r * .9) + "h" + f(1.8 * r) + "v" + f(1.8 * r) + "h" + f(-1.8 * r) + "z";
      case "diamond": { var q = r * 1.3; return "M" + f(x) + "," + f(y - q) + "L" + f(x + q) + "," + f(y) + "L" + f(x) + "," + f(y + q) + "L" + f(x - q) + "," + f(y) + "z"; }
      case "tri": return "M" + f(x) + "," + f(y - r * 1.25) + "L" + f(x + r * 1.2) + "," + f(y + r * .9) + "L" + f(x - r * 1.2) + "," + f(y + r * .9) + "z";
      case "trid": return "M" + f(x) + "," + f(y + r * 1.25) + "L" + f(x + r * 1.2) + "," + f(y - r * .9) + "L" + f(x - r * 1.2) + "," + f(y - r * .9) + "z";
      case "hex": { var s = ""; for (var i = 0; i < 6; i++) { var a = Math.PI / 3 * i; s += (i ? "L" : "M") + f(x + r * 1.15 * Math.cos(a)) + "," + f(y + r * 1.15 * Math.sin(a)); } return s + "z"; }
      case "plus": { var w = r * .42, a2 = r * 1.15; return "M" + f(x - w) + "," + f(y - a2) + "h" + f(2 * w) + "v" + f(a2 - w) + "h" + f(a2 - w) + "v" + f(2 * w) + "h" + f(-(a2 - w)) + "v" + f(a2 - w) + "h" + f(-2 * w) + "v" + f(-(a2 - w)) + "h" + f(-(a2 - w)) + "v" + f(-2 * w) + "h" + f(a2 - w) + "z"; }
    }
  }
  function starPath(x, y, r) {
    var s = "";
    for (var i = 0; i < 10; i++) { var a = -Math.PI / 2 + i * Math.PI / 5, rr = i % 2 ? r * .42 : r; s += (i ? "L" : "M") + (x + rr * Math.cos(a)).toFixed(1) + "," + (y + rr * Math.sin(a)).toFixed(1); }
    return s + "z";
  }
  function niceTicks(a, b, n) {
    var span = b - a, raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    var step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p, t = [];
    for (var v = Math.ceil(a / step - 1e-9) * step; v <= b + step * 1e-9; v += step) t.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toFixed(12));
    return { t: t, dec: Math.min(6, Math.max(0, -Math.floor(Math.log10(step) + 1e-9))) };
  }
  function baseExtent() {
    var b = { x0: 1e9, x1: -1e9, y0: 1e9, y1: -1e9 };
    pts.forEach(function (p) { var x = p[3 + px], y = p[3 + py]; if (x < b.x0) b.x0 = x; if (x > b.x1) b.x1 = x; if (y < b.y0) b.y0 = y; if (y > b.y1) b.y1 = y; });
    var dx = (b.x1 - b.x0) * .05, dy = (b.y1 - b.y0) * .05;
    b.x0 -= dx; b.x1 += dx; b.y0 -= dy; b.y1 += dy; return b;
  }
  function resetView() { base = baseExtent(); view = { x0: base.x0, x1: base.x1, y0: base.y0, y1: base.y1 }; }
  function ellipseAxes() { // returns "xy", "yx" or null for the current axes
    if (!P.ellipse) return null;
    if (px === 0 && py === 1) return "xy";
    if (px === 1 && py === 0) return "yx";
    return null;
  }

  var svg = $("pca"), wrap = $("plotwrap"), tip = $("tip");
  var M = { l: 60, r: 16, t: 12, b: 48 }, W = 800, H = 560, pw = 0, ph = 0, ringEl = null;
  function sx(v) { return M.l + (v - view.x0) / (view.x1 - view.x0) * pw; }
  function sy(v) { return M.t + ph - (v - view.y0) / (view.y1 - view.y0) * ph; }
  function ix(px_) { return view.x0 + (px_ - M.l) / pw * (view.x1 - view.x0); }
  function iy(py_) { return view.y0 + (M.t + ph - py_) / ph * (view.y1 - view.y0); }

  function drawPCA() {
    if (mode !== "pca") return;
    if (!view) resetView();
    W = Math.max(320, wrap.clientWidth); H = Math.round(Math.min(640, Math.max(400, W * .72)));
    pw = W - M.l - M.r; ph = H - M.t - M.b;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H); svg.setAttribute("height", H);
    var out = [];
    out.push('<defs><clipPath id="clip"><rect x="' + M.l + '" y="' + M.t + '" width="' + pw + '" height="' + ph + '"/></clipPath></defs>');
    out.push('<rect x="' + M.l + '" y="' + M.t + '" width="' + pw + '" height="' + ph + '" fill="var(--surface)" stroke="var(--line)"/>');
    var tx = niceTicks(view.x0, view.x1, Math.max(3, Math.floor(pw / 90))), ty = niceTicks(view.y0, view.y1, Math.max(3, Math.floor(ph / 70)));
    tx.t.forEach(function (v) { var X = sx(v); out.push('<line x1="' + X + '" x2="' + X + '" y1="' + M.t + '" y2="' + (M.t + ph) + '" stroke="' + (v === 0 ? "var(--line)" : "var(--grid)") + '"/><text class="tick" x="' + X + '" y="' + (M.t + ph + 16) + '" text-anchor="middle">' + v.toFixed(tx.dec) + '</text>'); });
    ty.t.forEach(function (v) { var Y = sy(v); out.push('<line y1="' + Y + '" y2="' + Y + '" x1="' + M.l + '" x2="' + (M.l + pw) + '" stroke="' + (v === 0 ? "var(--line)" : "var(--grid)") + '"/><text class="tick" x="' + (M.l - 8) + '" y="' + (Y + 3.5) + '" text-anchor="end">' + v.toFixed(ty.dec) + '</text>'); });
    out.push('<text class="axis-label" x="' + (M.l + pw / 2) + '" y="' + (H - 8) + '" text-anchor="middle">PC' + (px + 1) + ' · ' + P.pct[px] + '% of variance</text>');
    out.push('<text class="axis-label" transform="translate(16 ' + (M.t + ph / 2) + ') rotate(-90)" text-anchor="middle">PC' + (py + 1) + ' · ' + P.pct[py] + '% of variance</text>');
    out.push('<g clip-path="url(#clip)">');
    var paths = [];
    for (var c = 0; c < NC; c++) paths.push("");
    pts.forEach(function (p) {
      if (hidden[p[0]]) return;
      var X = sx(p[3 + px]), Y = sy(p[3 + py]);
      if (X < M.l - 6 || X > M.l + pw + 6 || Y < M.t - 6 || Y > M.t + ph + 6) return;
      paths[p[0]] += shapePath(SHAPES[p[0]], X, Y, 3.6);
    });
    for (c = 0; c < NC; c++) if (paths[c]) out.push('<path d="' + paths[c] + '" fill="var(--c' + c + ')" fill-opacity=".8"/>');
    var ea = ellipseAxes(), sxv = S.pcs[px], syv = S.pcs[py];
    if (showEll && ea) {
      var e = P.ellipse, poly = "";
      for (var i = 0; i <= 120; i++) {
        var th = i / 120 * 2 * Math.PI, ex = e.a * Math.cos(th), ey = e.b * Math.sin(th);
        var rx = Math.cos(e.ang) * ex - Math.sin(e.ang) * ey + e.cx, ry = Math.sin(e.ang) * ex + Math.cos(e.ang) * ey + e.cy;
        var X2 = ea === "xy" ? rx : ry, Y2 = ea === "xy" ? ry : rx;
        poly += (i ? "L" : "M") + sx(X2).toFixed(1) + "," + sy(Y2).toFixed(1);
      }
      out.push('<path d="' + poly + 'z" fill="var(--ink)" fill-opacity=".10" stroke="var(--ink)" stroke-width="1.6" stroke-dasharray="5 3"/>');
    }
    var SX = sx(sxv), SY = sy(syv);
    out.push('<path d="' + starPath(SX, SY, 10) + '" fill="var(--ink)" stroke="var(--surface)" stroke-width="2" stroke-linejoin="round"/>');
    out.push('<text x="' + (SX + 14) + '" y="' + (SY + 4) + '" font-size="13" font-weight="700" fill="var(--ink)" stroke="var(--surface)" stroke-width="3.5" paint-order="stroke" stroke-linejoin="round">__SAMPLE__</text>');
    out.push('<circle id="ring" r="9" fill="none" stroke="var(--ink)" stroke-width="2" visibility="hidden"/>');
    out.push('<rect id="band" fill="var(--ink)" fill-opacity=".08" stroke="var(--ink)" stroke-dasharray="4 3" visibility="hidden"/>');
    out.push('</g>');
    out.push('<rect id="hit" x="' + M.l + '" y="' + M.t + '" width="' + pw + '" height="' + ph + '" fill="transparent"/>');
    svg.innerHTML = out.join("");
    ringEl = $("ring");
    updateRail();
  }

  /* zoom + interaction */
  var MINSPAN = 0.0006;
  function setView(x0, x1, y0, y1) {
    if (x1 - x0 < MINSPAN) { var cx = (x0 + x1) / 2; x0 = cx - MINSPAN / 2; x1 = cx + MINSPAN / 2; }
    if (y1 - y0 < MINSPAN) { var cy = (y0 + y1) / 2; y0 = cy - MINSPAN / 2; y1 = cy + MINSPAN / 2; }
    view = { x0: x0, x1: x1, y0: y0, y1: y1 }; drawPCA();
  }
  function zoomAt(cx, cy, k) {
    setView(cx - (cx - view.x0) * k, cx + (view.x1 - cx) * k, cy - (cy - view.y0) * k, cy + (view.y1 - cy) * k);
  }
  $("zin").onclick = function () { zoomAt((view.x0 + view.x1) / 2, (view.y0 + view.y1) / 2, .5); };
  $("zout").onclick = function () { zoomAt((view.x0 + view.x1) / 2, (view.y0 + view.y1) / 2, 2); };
  $("zreset").onclick = function () { resetView(); drawPCA(); };
  $("zsample").onclick = function () {
    var cx = S.pcs[px], cy = S.pcs[py], hx, hy;
    if (P.ellipse) { hx = hy = P.ellipse.a * 6; } else { hx = (base.x1 - base.x0) * .03; hy = (base.y1 - base.y0) * .03; }
    setView(cx - hx, cx + hx, cy - hy, cy + hy);
  };

  var drag = null;
  function local(e) { var r = svg.getBoundingClientRect(); return { x: (e.clientX - r.left) * (W / r.width), y: (e.clientY - r.top) * (H / r.height) }; }
  svg.addEventListener("pointerdown", function (e) {
    if (e.target.id !== "hit") return;
    var p = local(e); drag = { x: p.x, y: p.y, moved: false, pan: e.shiftKey, v: { x0: view.x0, x1: view.x1, y0: view.y0, y1: view.y1 } };
    svg.setPointerCapture(e.pointerId); tip.hidden = true;
  });
  svg.addEventListener("pointermove", function (e) {
    var p = local(e);
    if (drag) {
      var dx = p.x - drag.x, dy = p.y - drag.y;
      if (Math.abs(dx) + Math.abs(dy) > 4) drag.moved = true;
      if (drag.moved && drag.pan) {
        var ddx = -dx / pw * (drag.v.x1 - drag.v.x0), ddy = dy / ph * (drag.v.y1 - drag.v.y0);
        view = { x0: drag.v.x0 + ddx, x1: drag.v.x1 + ddx, y0: drag.v.y0 + ddy, y1: drag.v.y1 + ddy };
        if (!drag.raf) drag.raf = requestAnimationFrame(function () { drag && (drag.raf = 0); drawPCA(); });
      } else if (drag.moved) {
        var b = $("band"); if (b) {
          var x0 = Math.max(M.l, Math.min(drag.x, p.x)), x1 = Math.min(M.l + pw, Math.max(drag.x, p.x)), y0 = Math.max(M.t, Math.min(drag.y, p.y)), y1 = Math.min(M.t + ph, Math.max(drag.y, p.y));
          b.setAttribute("x", x0); b.setAttribute("y", y0); b.setAttribute("width", x1 - x0); b.setAttribute("height", y1 - y0); b.setAttribute("visibility", "visible");
        }
      }
      return;
    }
    hover(p, e);
  });
  function endDrag(e) {
    if (!drag) return;
    var p = local(e), d = drag; drag = null;
    if (!d.moved) { var k = e.shiftKey ? 2 : .5; zoomAt(ix(p.x), iy(p.y), k); return; }
    if (d.pan) return;
    var xa = Math.max(M.l, Math.min(d.x, p.x)), xb = Math.min(M.l + pw, Math.max(d.x, p.x)), ya = Math.max(M.t, Math.min(d.y, p.y)), yb = Math.min(M.t + ph, Math.max(d.y, p.y));
    if (xb - xa < 8 || yb - ya < 8) { drawPCA(); return; }
    setView(ix(xa), ix(xb), iy(yb), iy(ya));
  }
  svg.addEventListener("pointerup", endDrag);
  svg.addEventListener("pointercancel", function () { drag = null; drawPCA(); });
  svg.addEventListener("pointerleave", function () { if (!drag) { tip.hidden = true; if (ringEl) ringEl.setAttribute("visibility", "hidden"); } });

  function hover(p, e) {
    if (p.x < M.l || p.x > M.l + pw || p.y < M.t || p.y > M.t + ph) { tip.hidden = true; return; }
    var best = null, bd = 15 * 15;
    var sxp = sx(S.pcs[px]), syp = sy(S.pcs[py]), ds = (sxp - p.x) * (sxp - p.x) + (syp - p.y) * (syp - p.y);
    if (ds < 16 * 16) best = { me: true, x: sxp, y: syp };
    else pts.forEach(function (q) {
      if (hidden[q[0]]) return;
      var X = sx(q[3 + px]), Y = sy(q[3 + py]), d = (X - p.x) * (X - p.x) + (Y - p.y) * (Y - p.y);
      if (d < bd) { bd = d; best = { q: q, x: X, y: Y }; }
    });
    if (!best) { tip.hidden = true; ringEl.setAttribute("visibility", "hidden"); return; }
    ringEl.setAttribute("cx", best.x); ringEl.setAttribute("cy", best.y); ringEl.setAttribute("r", best.me ? 15 : 8); ringEl.setAttribute("visibility", "visible");
    tip.textContent = "";
    var head = el("b", null, best.me ? "__SAMPLE__" : P.subs[best.q[1]]);
    var sub = el("span", null, best.me ? "Case sample" : CONT[best.q[0]].label + " · " + best.q[2]);
    var v = best.me ? S.pcs : best.q.slice(3);
    var line = el("div", "m", "PC" + (px + 1) + " " + v[px].toFixed(4) + "  PC" + (py + 1) + " " + v[py].toFixed(4));
    tip.appendChild(head); tip.appendChild(sub); tip.appendChild(line); tip.hidden = false;
    var tw = tip.offsetWidth, th = tip.offsetHeight;
    var r = svg.getBoundingClientRect(), scale = r.width / W;
    var left = best.x * scale + 14, top = best.y * scale - th - 12;
    if (left + tw > r.width) left = best.x * scale - tw - 14;
    if (top < 0) top = best.y * scale + 16;
    tip.style.left = left + "px"; tip.style.top = top + "px";
  }

  /* rail controls */
  function buildPCButtons(host, which) {
    host.textContent = "";
    for (var i = 0; i < NPC; i++) (function (i) {
      var b = el("button", "chip", "PC" + (i + 1)); b.type = "button";
      b.setAttribute("aria-pressed", (which === "x" ? px : py) === i);
      b.onclick = function () {
        if (which === "x") { if (i === py) py = px; px = i; } else { if (i === px) px = py; py = i; }
        buildPCButtons($("xpcs"), "x"); buildPCButtons($("ypcs"), "y"); resetView(); drawPCA();
      };
      host.appendChild(b);
    })(i);
  }
  $("swap").onclick = function () { var t = px; px = py; py = t; buildPCButtons($("xpcs"), "x"); buildPCButtons($("ypcs"), "y"); resetView(); drawPCA(); };
  $("ell").onchange = function () { showEll = this.checked; drawPCA(); };

  function swatch(c, sz) {
    var s = document.createElementNS("http://www.w3.org/2000/svg", "svg"); s.setAttribute("width", sz); s.setAttribute("height", sz); s.setAttribute("viewBox", "-8 -8 16 16");
    var p = document.createElementNS("http://www.w3.org/2000/svg", "path"); p.setAttribute("d", shapePath(SHAPES[c], 0, 0, 5.6)); p.setAttribute("fill", "var(--c" + c + ")"); s.appendChild(p); return s;
  }
  function buildLegend() {
    var host = $("legend"); host.textContent = "";
    CONT.forEach(function (c, i) {
      var b = el("button"); b.type = "button"; b.setAttribute("aria-pressed", !hidden[i]);
      b.appendChild(swatch(i, 16)); b.appendChild(el("span", null, c.label)); b.appendChild(el("span", "n", String(c.n)));
      b.onclick = function () { hidden[i] = !hidden[i]; buildLegend(); drawPCA(); };
      host.appendChild(b);
    });
  }
  $("showall").onclick = function () { hidden = {}; buildLegend(); drawPCA(); };

  function updateRail() {
    var ea = ellipseAxes(), cb = $("ell");
    cb.disabled = !ea;
    $("ellhint").textContent = !P.ellipse ? "No confidence ellipse was found for this run." : ea ? "Around __SAMPLE__, from smartpca (ellipse run, PC1 and PC2)." : "Available for PC1 vs PC2 only. Switch the axes to see it.";
    var ro = $("readout"); ro.textContent = "";
    [["PC" + (px + 1), S.pcs[px].toFixed(4)], ["PC" + (py + 1), S.pcs[py].toFixed(4)]].forEach(function (r) { ro.appendChild(el("dt", null, r[0])); ro.appendChild(el("dd", null, r[1])); });
    // nearest populations by centroid distance in the two displayed PCs
    var acc = {};
    pts.forEach(function (q) { var a = acc[q[1]] || (acc[q[1]] = { n: 0, x: 0, y: 0, c: q[0] }); a.n++; a.x += q[3 + px]; a.y += q[3 + py]; });
    var list = Object.keys(acc).map(function (k) { var a = acc[k], dx = a.x / a.n - S.pcs[px], dy = a.y / a.n - S.pcs[py]; return { name: P.subs[k], c: a.c, d: Math.sqrt(dx * dx + dy * dy) }; }).sort(function (a, b) { return a.d - b.d; }).slice(0, 3);
    var ul = $("near"); ul.textContent = "";
    var head = el("li"); head.style.color = "var(--ink-2)"; head.textContent = "Nearest populations in these two PCs"; ul.appendChild(head);
    list.forEach(function (r) { var li = el("li"); li.appendChild(swatch(r.c, 14)); li.appendChild(el("span", null, r.name)); li.appendChild(el("span", "d", r.d.toFixed(4))); ul.appendChild(li); });
  }

  /* ================= ADMIXTURE ================= */
  var A = D.admix || { byK: {}, cv: {}, pops: [] }, Ks = Object.keys(A.byK).map(Number).sort(function (a, b) { return a - b; });
  var K = Ks[0], bestK = Ks[0];
  Ks.forEach(function (k) { if (A.cv[k] < A.cv[bestK]) bestK = k; });
  K = bestK;
  var sel = {}; A.pops.forEach(function (p, i) { sel[i] = true; });
  var admView = "bars", query = "";

  function components(k) {
    var avg = A.byK[k].avg, cm = [];
    for (var c = 0; c < NC; c++) { cm.push({ w: 0, s: new Array(k).fill(0) }); }
    A.pops.forEach(function (p, i) { cm[p.cont].w += p.n; for (var j = 0; j < k; j++) cm[p.cont].s[j] += avg[i][j] * p.n; });
    var dom = [];
    for (var j = 0; j < k; j++) {
      var bc = 0, bv = -1;
      for (c = 0; c < NC; c++) { var v = cm[c].s[j] / cm[c].w; if (v > bv) { bv = v; bc = c; } }
      dom.push({ j: j, cont: bc, str: bv });
    }
    var groups = {}; dom.forEach(function (d) { (groups[d.cont] = groups[d.cont] || []).push(d); });
    var out = [];
    for (c = 0; c < NC; c++) {
      var g = (groups[c] || []).sort(function (a, b) { return b.str - a.str; });
      g.forEach(function (d, r) { out.push({ j: d.j, cont: c, rank: r, label: CONT[c].label + (g.length > 1 ? " " + (r + 1) : ""), css: compCss(c, r) }); });
    }
    return out;
  }
  function compCss(c, r) {
    var base = "var(--c" + c + ")", tint = "color-mix(in srgb, var(--c" + c + ") 40%, var(--surface))";
    switch (r % 4) {
      case 0: return base;
      case 1: return "repeating-linear-gradient(135deg," + base + " 0 3px," + tint + " 3px 6px)";
      case 2: return tint;
      default: return "repeating-linear-gradient(45deg," + tint + " 0 3px," + base + " 3px 6px)";
    }
  }

  function buildKs() {
    var host = $("ks"); host.textContent = "";
    Ks.forEach(function (k) {
      var b = el("button", "kbtn"); b.type = "button"; b.setAttribute("aria-pressed", k === K);
      b.appendChild(el("b", null, "K=" + k)); b.appendChild(el("small", null, "CV " + A.cv[k].toFixed(4)));
      if (k === bestK) b.appendChild(el("span", "low", "lowest CV"));
      b.onclick = function () { K = k; buildKs(); drawAdm(); };
      host.appendChild(b);
    });
  }
  function buildPicker() {
    var host = $("plist"); host.textContent = "";
    CONT.forEach(function (c, ci) {
      var idx = []; A.pops.forEach(function (p, i) { if (p.cont === ci && p.name.toLowerCase().indexOf(query) > -1) idx.push(i); });
      if (!idx.length) return;
      var g = el("div", "pgroup"), gh = el("div", "gh");
      var all = el("input"); all.type = "checkbox"; all.id = "g" + ci; all.setAttribute("aria-label", "Toggle " + c.label);
      var on = idx.filter(function (i) { return sel[i]; }).length; all.checked = on === idx.length; all.indeterminate = on > 0 && on < idx.length;
      all.onchange = function () { idx.forEach(function (i) { sel[i] = all.checked; }); buildPicker(); drawAdm(); };
      var sw = el("span", "sw"); sw.style.background = "var(--c" + ci + ")";
      gh.appendChild(all); gh.appendChild(sw); gh.appendChild(el("span", null, c.label)); g.appendChild(gh);
      idx.forEach(function (i) {
        var lab = el("label", "p"), cb = el("input"); cb.type = "checkbox"; cb.checked = !!sel[i];
        cb.onchange = function () { sel[i] = cb.checked; buildPicker(); drawAdm(); };
        lab.appendChild(cb); lab.appendChild(el("span", null, A.pops[i].name)); lab.appendChild(el("span", "n", String(A.pops[i].n))); g.appendChild(lab);
      });
      host.appendChild(g);
    });
    if (!host.firstChild) host.appendChild(el("p", "hint", "No population matches “" + query + "”."));
  }
  $("psearch").oninput = function () { query = this.value.trim().toLowerCase(); buildPicker(); };
  $("pall").onclick = function () { A.pops.forEach(function (p, i) { sel[i] = true; }); buildPicker(); drawAdm(); };
  $("pnone").onclick = function () { A.pops.forEach(function (p, i) { sel[i] = false; }); buildPicker(); drawAdm(); };
  $("v-bars").onclick = function () { admView = "bars"; $("v-bars").setAttribute("aria-pressed", true); $("v-table").setAttribute("aria-pressed", false); drawAdm(); };
  $("v-table").onclick = function () { admView = "table"; $("v-table").setAttribute("aria-pressed", true); $("v-bars").setAttribute("aria-pressed", false); drawAdm(); };

  var tip2 = $("tip2");
  function showTip(e, title, sub, rows, hi) {
    tip2.textContent = "";
    tip2.appendChild(el("b", null, title)); if (sub) tip2.appendChild(el("span", null, sub));
    rows.forEach(function (r) {
      var d = el("div", "k"), i = el("i"); i.style.background = r.css; d.appendChild(i);
      var t = el("span", "m", (r.j === hi ? "▸ " : "") + pct(r.v) + "  " + r.label); if (r.j === hi) t.style.fontWeight = "700"; d.appendChild(t); tip2.appendChild(d);
    });
    tip2.hidden = false;
    var host = tip2.parentElement.getBoundingClientRect(), tw = tip2.offsetWidth, th = tip2.offsetHeight;
    var x = e.clientX - host.left + 14, y = e.clientY - host.top + 14;
    if (x + tw > host.width) x = e.clientX - host.left - tw - 14;
    if (y + th > host.height) y = e.clientY - host.top - th - 14;
    tip2.style.left = Math.max(0, x) + "px"; tip2.style.top = Math.max(0, y) + "px";
  }
  function rowTip(e, name, sub, q, comps, hi) {
    var rows = comps.map(function (c) { return { j: c.j, v: q[c.j], label: c.label, css: c.css }; }).filter(function (r) { return r.v >= .01 || r.j === hi; }).sort(function (a, b) { return b.v - a.v; });
    showTip(e, name, sub, rows, hi);
  }

  function drawAdm() {
    if (mode !== "adm") return;
    var comps = components(K), sample = A.byK[K].sample, avg = A.byK[K].avg;
    // summary of the case sample, grouped by the continent each component belongs to
    var byC = {}; comps.forEach(function (c) { byC[c.cont] = (byC[c.cont] || 0) + sample[c.j]; });
    var sm = $("summary"); sm.textContent = ""; sm.appendChild(el("h3", null, "__SAMPLE__ at K=" + K));
    Object.keys(byC).map(function (c) { return { c: +c, v: byC[c] }; }).sort(function (a, b) { return b.v - a.v; }).filter(function (r) { return r.v >= .02; }).forEach(function (r) {
      var s = el("span", "part"), i = el("i"); i.style.background = "var(--c" + r.c + ")"; s.appendChild(i); s.appendChild(el("span", "mono", pct(r.v))); s.appendChild(document.createTextNode(CONT[r.c].label)); sm.appendChild(s);
    });
    var cl = $("comps"); cl.textContent = "";
    comps.forEach(function (c) { var s = el("span"), i = el("i"); i.style.background = c.css; s.appendChild(i); s.appendChild(document.createTextNode(c.label)); cl.appendChild(s); });

    var body = $("admbody"); body.textContent = ""; tip2.hidden = true;
    var chosen = []; A.pops.forEach(function (p, i) { if (sel[i]) chosen.push(i); });
    if (admView === "bars") {
      var box = el("div", "bars");
      function bar(q, big) {
        var b = el("div", "bar");
        comps.forEach(function (c) {
          var v = q[c.j]; if (v < .002) return;
          var seg = el("i"); seg.style.width = (v * 100).toFixed(3) + "%"; seg.style.background = c.css;
          seg.setAttribute("aria-label", c.label + " " + pct(v));
          b.appendChild(seg);
        });
        return b;
      }
      var me = el("div", "brow me"); me.appendChild(el("div", "lab", "__SAMPLE__")); var mb = bar(sample, true); me.appendChild(mb); box.appendChild(me);
      mb.addEventListener("pointermove", function (e) { var t = e.target; var hi = null; showTip(e, "__SAMPLE__", "Case sample · K=" + K, comps.map(function (c) { return { j: c.j, v: sample[c.j], label: c.label, css: c.css }; }).filter(function (r) { return r.v >= .01; }).sort(function (a, b) { return b.v - a.v; }), hi); });
      mb.addEventListener("pointerleave", function () { tip2.hidden = true; });
      if (!chosen.length) box.appendChild(el("p", "empty", "No reference populations selected. Choose some on the left."));
      var lastC = -1;
      chosen.forEach(function (i) {
        var p = A.pops[i];
        if (p.cont !== lastC) { box.appendChild(el("div", "grp", CONT[p.cont].label)); lastC = p.cont; }
        var r = el("div", "brow"), lab = el("div", "lab", p.name); lab.title = p.name + " (" + p.n + " individuals)"; lab.appendChild(el("small", null, String(p.n)));
        var b = bar(avg[i]); r.appendChild(lab); r.appendChild(b); box.appendChild(r);
        b.addEventListener("pointermove", function (e) {
          var idx = Array.prototype.indexOf.call(b.children, e.target), hi = null;
          if (idx > -1) { var shown = comps.filter(function (c) { return avg[i][c.j] >= .002; }); hi = shown[idx] ? shown[idx].j : null; }
          rowTip(e, p.name, CONT[p.cont].label + " · mean of " + p.n + " individuals", avg[i], comps, hi);
        });
        b.addEventListener("pointerleave", function () { tip2.hidden = true; });
      });
      body.appendChild(box);
    } else {
      var wrapT = el("div", "tablewrap"), t = el("table"), thead = el("thead"), hr = el("tr");
      var th0 = el("th", null, "Population"); hr.appendChild(th0); hr.appendChild(el("th", null, "n"));
      comps.forEach(function (c) { var th = el("th"), i = el("i"); i.style.background = c.css; th.appendChild(i); th.appendChild(document.createTextNode(c.label)); hr.appendChild(th); });
      thead.appendChild(hr); t.appendChild(thead);
      var tb = el("tbody");
      function trow(name, n, q, cls) {
        var tr = el("tr", cls); tr.appendChild(el("th", null, name)); tr.appendChild(el("td", "dim", n));
        comps.forEach(function (c) { var td = el("td", q[c.j] < .005 ? "dim" : null, (q[c.j] * 100).toFixed(1)); tr.appendChild(td); }); tb.appendChild(tr);
      }
      trow("__SAMPLE__", "1", sample, "me");
      chosen.forEach(function (i) { trow(A.pops[i].name, String(A.pops[i].n), avg[i]); });
      t.appendChild(tb); wrapT.appendChild(t); body.appendChild(wrapT);
      body.appendChild(el("p", "hint", "Values are percent ancestry. Population rows are means over the reference individuals."));
    }
    var cvs = Ks.map(function (k) { return A.cv[k]; }), spread = Math.max.apply(null, cvs) - Math.min.apply(null, cvs);
    $("admnote").textContent = "Each bar is the mean ancestry of a reference population, except __SAMPLE__, which is one sample" + (m.coverage != null ? " at " + m.coverage.toFixed(2) + "× coverage" : "") + ", so its proportions are indicative. " +
      "Cross-validation error varies by only " + spread.toFixed(3) + " across K=" + Ks[0] + " to " + Ks[Ks.length - 1] + ", so no K is clearly best. " +
      "Component colours follow the continent where each component is highest; components that share a continent are told apart by pattern.";
  }

  /* ---------- boot ---------- */
  buildPCButtons($("xpcs"), "x"); buildPCButtons($("ypcs"), "y"); buildLegend(); buildKs(); buildPicker();
  resetView(); drawPCA();
  if (window.ResizeObserver) new ResizeObserver(function () { if (mode === "pca") drawPCA(); }).observe(wrap); else window.addEventListener("resize", function () { drawPCA(); });
})();
</script>
'''

if __name__ == "__main__":
    main()
