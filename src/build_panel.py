#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = []
# ///
"""Step 5 - build the a-priori gene panel from the patient's HPO terms.

Ranks DISEASES against the patient's phenotype, then takes their genes. Ranking genes
directly does not work: four general dysmorphology terms union to about half of every gene
HPO knows, which is not a panel. Disease-level ranking weights the *combination* of terms.

Runs before any variant is seen and writes its output with a timestamp, so the claim that
the phenotype ranking was not tuned against the answer is checkable rather than asserted.
Needs no VCF and no annotation - run it while the annotation is still going.

    ./src/build_panel.py --terms config/patient_hpo.txt

Outputs (default <repo>/results):
    panel.tsv        ranked genes with scores, supporting disease and coordinates
    panel.meta.json  full record INCLUDING the patient's HPO terms - patient material

Neither the inputs nor panel.meta.json may be committed: they carry the proband's clinical
features. The challenge's participant reference document requires that nothing capable of
re-identifying the patient or family be published beyond what the family has already made
public themselves, so this script checks on every run that those files are untracked and
refuses to proceed if any of them is in the index.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Root of "Phenotypic abnormality". Terms outside this subtree (inheritance, clinical
#: course, modifiers) are not phenotype and would distort the information content.
PHENOTYPIC_ABNORMALITY = "HP:0000118"

#: Diseases scoring below this fraction of the best disease contribute no genes. Not a
#: significance threshold - there is no null model here - just a cut that keeps the panel
#: from trailing off into thousands of genes sharing one grandparent term.
MIN_REL_SCORE = 0.40


def log(msg):
    print(msg, flush=True)


#: Files that hold, or can reconstruct, the proband's phenotype. Kept out of version control
#: as a matter of the challenge's data-privacy terms, and checked rather than assumed: a
#: .gitignore entry stops a NEW file being added but does nothing about one already committed,
#: and `git add -f` overrides it without complaint.
PRIVATE_PHENOTYPE_FILES = ("config/patient_hpo.txt",
                           "config/family_history.txt",
                           "results/panel.meta.json")


def check_phenotype_privacy(root=ROOT):
    """Refuse to run if a phenotype-bearing file is tracked by git.

    The challenge's participant reference document states that the clinical information is
    real, from a living patient and family who consented to share it, and that nothing which
    could re-identify them beyond what the family has already published may appear in any
    output of this work. Asking git what it tracks is the only reliable way to check that.
    """
    if not (root / ".git").is_dir():
        log("  phenotype privacy: not a git repository, check skipped")
        return
    tracked = [rel for rel in PRIVATE_PHENOTYPE_FILES
               if subprocess.run(["git", "-C", str(root), "ls-files", "--error-unmatch", rel],
                                 capture_output=True, text=True).returncode == 0]
    if tracked:
        die("phenotype file(s) tracked by git: " + ", ".join(tracked)
            + "\n  These carry the proband's clinical features, which the challenge's data"
              " terms keep out of anything published."
              "\n  Take them out of the index before continuing:"
            + "".join(f"\n    git rm --cached {t}" for t in tracked))
    log(f"  phenotype privacy: {len(PRIVATE_PHENOTYPE_FILES)} file(s) checked, none tracked")


def die(msg):
    sys.exit(f"ERROR: {msg}")


# --------------------------------------------------------------------------- hp.obo

def parse_obo(path):
    """Return (terms, alt_to_id, version). terms: id -> name/parents/obsolete/replaced_by."""
    terms, alt = {}, {}
    cur = None
    version = "unknown"
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.rstrip()
        if line.startswith("data-version:"):
            version = line.split(":", 1)[1].strip()
        elif line == "[Term]":
            cur = {"id": None, "name": "", "parents": set(),
                   "obsolete": False, "replaced_by": None}
        elif line.startswith("["):          # [Typedef] and friends end the term block
            cur = None
        elif cur is not None and line:
            key, _, val = line.partition(": ")
            if key == "id" and val.startswith("HP:"):
                cur["id"] = val
                terms[val] = cur
            elif key == "name":
                cur["name"] = val
            elif key == "is_a":
                cur["parents"].add(val.split(" ! ")[0].strip())
            elif key == "alt_id" and cur["id"]:
                alt[val] = cur["id"]
            elif key == "is_obsolete" and val == "true":
                cur["obsolete"] = True
            elif key == "replaced_by":
                cur["replaced_by"] = val.split(" ! ")[0].strip()
    if not terms:
        die(f"no [Term] stanzas in {path} - is this really hp.obo?")
    return terms, alt, version


def make_ancestor_fn(terms):
    """Ancestors of a term, inclusive of itself. Memoised; the DAG is walked once."""
    @lru_cache(maxsize=None)
    def anc(tid):
        out = {tid}
        for p in terms.get(tid, {}).get("parents", ()):   # unknown ids contribute nothing
            out |= anc(p)
        return frozenset(out)
    return anc


def canonicalise(raw, terms, alt):
    """Map an id through alt_id and obsolete -> replaced_by. Returns (id|None, note)."""
    tid = alt.get(raw, raw)
    note = "alt_id" if tid != raw else ""
    t = terms.get(tid)
    if t is None:
        return None, "unknown id"
    if t["obsolete"]:
        if t["replaced_by"] and t["replaced_by"] in terms:
            return t["replaced_by"], "obsolete -> replaced_by"
        return None, "obsolete, no replacement"
    return tid, note


# ------------------------------------------------------------------- phenotype.hpoa

def load_hpoa(path, anc, keep_subtree):
    """disease_id -> {name, terms(propagated)}. Keeps aspect P, drops NOT-qualified rows."""
    diseases, version, cols = {}, "unknown", None
    dropped_not = dropped_aspect = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if line.lower().startswith("#version"):
                    version = line.split(":", 1)[1].strip()
                continue
            f = line.split("\t")
            if cols is None:
                if f[0] != "database_id":
                    die(f"{path}: expected a 'database_id' header line, got {f[0]!r}")
                cols = {name: i for i, name in enumerate(f)}
                for need in ("database_id", "disease_name", "hpo_id", "aspect"):
                    if need not in cols:
                        die(f"{path}: header has no {need!r} column")
                continue
            if len(f) <= max(cols.values()):
                continue
            if f[cols["aspect"]] != "P":
                dropped_aspect += 1
                continue
            # A NOT row says the disease does NOT have this feature. Counting it would
            # score the patient against the negation of the disease.
            if "qualifier" in cols and f[cols["qualifier"]].strip().upper() == "NOT":
                dropped_not += 1
                continue
            tid = f[cols["hpo_id"]].strip()
            if not tid.startswith("HP:"):
                continue
            d = diseases.setdefault(f[cols["database_id"]],
                                    {"name": f[cols["disease_name"]], "terms": set()})
            d["terms"] |= set(anc(tid)) & keep_subtree
    if cols is None:
        die(f"{path}: no header line found")
    diseases = {k: {"name": v["name"], "terms": frozenset(v["terms"])}
                for k, v in diseases.items() if v["terms"]}
    if not diseases:
        die(f"{path}: no usable disease annotations after filtering")
    log(f"  diseases: {len(diseases):,} (dropped {dropped_aspect:,} non-P rows, "
        f"{dropped_not:,} NOT rows)")
    return diseases, version


def information_content(diseases):
    """IC(t) = -log(P(t)), P estimated over DISEASES, on already-propagated term sets."""
    n = len(diseases)
    freq = defaultdict(int)
    for d in diseases.values():
        for t in d["terms"]:
            freq[t] += 1
    return {t: -math.log(c / n) for t, c in freq.items()}


# ----------------------------------------------------------------------- similarity

def build_sim_rows(patient, all_terms, ic, anc):
    """sim[p][t] = Resnik similarity (IC of the most informative common ancestor).

    Precomputed per patient term over every term HPO uses, so scoring 13k diseases is
    dictionary lookups rather than tens of millions of set intersections.
    """
    rows = {}
    for p in patient:
        ap = anc(p)
        row = {}
        for t in all_terms:
            common = ap & anc(t)
            if common:
                v = max(ic.get(c, 0.0) for c in common)
                if v > 0:
                    row[t] = v
        rows[p] = row
    return rows


def score_disease(patient, disease_terms, rows, best_for_term):
    """Best-match averages in both directions.

    patient->disease asks 'does this disease explain what we see'. The reverse direction
    penalises a disease for features the child was not reported to have, which with an
    incomplete term list punishes exactly the well-described syndromes. Ranking uses the
    forward score; the symmetric average is reported alongside it as a tiebreak.
    """
    if not patient or not disease_terms:
        return 0.0, 0.0
    fwd = sum(max((rows[p].get(d, 0.0) for d in disease_terms), default=0.0)
              for p in patient) / len(patient)
    rev = sum(best_for_term.get(d, 0.0) for d in disease_terms) / len(disease_terms)
    return fwd, (fwd + rev) / 2


# ---------------------------------------------------------------------------- genes

def load_disease_genes(path):
    """disease_id -> {gene symbols}, from genes_to_phenotype.txt."""
    mapping = defaultdict(set)
    with path.open(encoding="utf-8", errors="replace") as fh:
        header = next(fh).rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        gi = cols.get("gene_symbol", 1)
        di = cols.get("disease_id", 5)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(gi, di) and f[di]:
                mapping[f[di]].add(f[gi])
    if not mapping:
        die(f"{path}: no disease->gene rows parsed")
    return mapping


def load_coords(path):
    """gene -> (chrom, start, end, strand) from refs/gene_coords.tsv."""
    if not path.exists():
        log(f"  ! {path} missing - panel will have no coordinates")
        return {}
    out = {}
    for i, line in enumerate(path.read_text().splitlines()):
        f = line.rstrip("\n").split("\t")
        if i == 0 and f[0] == "gene":
            continue
        if len(f) >= 5:
            out[f[0]] = (f[1], f[2], f[3], f[4])
    return out


def load_checklist(mva_tsv, extra_tsv):
    """MVA genes that appear whatever the phenotype ranking says."""
    genes = set()
    if mva_tsv.exists():
        for i, line in enumerate(mva_tsv.read_text().splitlines()):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if i == 0 and f[0] == "gene":
                continue
            genes.add(f[0].strip())
    else:
        log(f"  ! {mva_tsv} missing - checklist comes only from {extra_tsv.name}")
    if extra_tsv.exists():
        for line in extra_tsv.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if f[0].strip() != "gene":
                genes.add(f[0].strip())
    return genes


# ----------------------------------------------------------------------------- main

HP_ID = re.compile(r"HP:\d{7}")


def read_terms(path):
    """One term per line, in whatever column order the clinical table used.

    The id is found ANYWHERE on the line, so a row pasted straight out of the challenge's
    clinical-features table works in either order ("HP:0000252\tMicrocephaly" or
    "Microcephaly, HP:0000252"). Everything else on the line is kept as a label and
    cross-checked against the ontology later: a correct name next to a mistyped id is the
    transcription error that silently rebuilds the panel around the wrong phenotype.
    '#' comments out the rest of a line.
    """
    if not path.exists():
        die(f"{path} not found. Create it with one HP:####### per line "
            f"(any column order, text around the id is kept as a label).")
    out, bad = [], []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        found = HP_ID.findall(line)
        if not found:
            bad.append(f"line {n}: {line[:60]!r}")
            continue
        if len(found) > 1:
            die(f"{path} line {n} has {len(found)} HPO ids on it ({', '.join(found)}) - "
                f"put one term per line")
        label = HP_ID.sub("", line).strip(" \t,;:|-").strip()
        out.append((found[0], label))
    if bad:
        # Loud: a row the parser skipped is a feature the child has and the panel does not.
        die(f"{path}: no HPO id on {len(bad)} non-empty line(s):\n    "
            + "\n    ".join(bad)
            + "\n  Comment them out with '#' if they are headers or free text.")
    if not out:
        die(f"{path} has no terms in it")
    return out


def label_disagrees(label, official):
    """True when a supplied label shares no word stem with the ontology's name.

    Deliberately weak: synonyms are normal and this must not cry wolf. It fires on the
    case that matters - a plausible name sitting next to an id for something unrelated.
    """
    if not label:
        return False
    def words(x):
        return [w for w in re.split(r"[^a-z0-9]+", x.lower()) if len(w) >= 4]
    a, b = words(label), words(official)
    if not a or not b:
        return False
    return not any(x.startswith(y[:4]) or y.startswith(x[:4]) for x in a for y in b)


def write_atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--terms", default=ROOT / "config" / "patient_hpo.txt", type=Path,
                    help="patient HPO terms, one HP:####### per line")
    ap.add_argument("--hpo-dir", default=ROOT / "refs" / "hpo", type=Path)
    ap.add_argument("--gene-coords", default=ROOT / "refs" / "gene_coords.tsv", type=Path)
    ap.add_argument("--mva-genes", default=ROOT / "refs" / "mva_genes.tsv", type=Path)
    ap.add_argument("--extra-genes", default=ROOT / "config" / "mva_genes_extra.tsv",
                    type=Path)
    ap.add_argument("--outdir", default=ROOT / "results", type=Path)
    ap.add_argument("--leave-one-out", action="store_true",
                    help="re-rank with each patient term removed in turn and report how "
                         "the top of the list moves. A highly specific term (a rare, "
                         "deeply nested feature) carries most of the information content "
                         "and can carry the ranking on its own - this says whether it "
                         "did, which matters when that feature might be incidental "
                         "rather than part of the syndrome")
    ap.add_argument("--min-rel-score", type=float, default=MIN_REL_SCORE,
                    help=f"disease cut as a fraction of the top score "
                         f"(default {MIN_REL_SCORE})")
    args = ap.parse_args()

    # Before anything is read or written: the phenotype files must not be in version
    # control. Enforced here rather than trusted to .gitignore - see the module docstring.
    check_phenotype_privacy()

    obo = args.hpo_dir / "hp.obo"
    hpoa = args.hpo_dir / "phenotype.hpoa"
    g2p = args.hpo_dir / "genes_to_phenotype.txt"
    for p in (obo, hpoa, g2p):
        if not p.exists():
            die(f"{p} not found - run ./src/download_refs.py first")

    log("[1/4] ontology")
    terms, alt, hpo_version = parse_obo(obo)
    anc = make_ancestor_fn(terms)
    keep = {t for t in terms if PHENOTYPIC_ABNORMALITY in anc(t)}
    log(f"  hp.obo {hpo_version}: {len(terms):,} terms, "
        f"{len(keep):,} under Phenotypic abnormality")

    patient, notes, unmapped, mismatched = [], {}, [], []
    for raw, label in read_terms(args.terms):
        tid, note = canonicalise(raw, terms, alt)
        if tid is None or tid not in keep:
            unmapped.append(f"{raw} ({note or 'outside Phenotypic abnormality'})")
            continue
        if note:
            notes[raw] = f"{note} -> {tid}"
        if label_disagrees(label, terms[tid]["name"]):
            mismatched.append(f"{raw} written as {label!r}, HPO calls it "
                              f"{terms[tid]['name']!r}")
        patient.append(tid)
    patient = sorted(set(patient))
    if unmapped:
        # Loud, not silent: a dropped term shortens the phenotype and the panel still
        # looks perfectly plausible.
        log(f"  ! {len(unmapped)} term(s) not usable: {', '.join(unmapped)}")
    if mismatched:
        # Advisory, not fatal - synonyms are normal. But a right-looking name beside a
        # wrong id rebuilds the panel around a phenotype the child does not have.
        log(f"  ! {len(mismatched)} label(s) disagree with the ontology - check for a "
            f"mistyped id:")
        for m in mismatched:
            log(f"      {m}")
    if not patient:
        die("none of the patient's terms could be resolved")
    log(f"  patient terms in use: {len(patient)}")
    for t in patient:
        log(f"    {t}  {terms[t]['name']}")

    log("[2/4] disease annotations")
    diseases, hpoa_version = load_hpoa(hpoa, anc, keep)
    ic = information_content(diseases)

    log("[3/4] ranking")
    all_terms = set().union(*(d["terms"] for d in diseases.values()))
    rows_sim = build_sim_rows(patient, all_terms, ic, anc)
    best_for_term = {t: max(rows_sim[p].get(t, 0.0) for p in patient) for t in all_terms}
    scored = []
    for did, d in diseases.items():
        fwd, sym = score_disease(patient, d["terms"], rows_sim, best_for_term)
        if fwd > 0:
            scored.append((fwd, sym, did, d["name"]))
    scored.sort(reverse=True)
    if not scored:
        die(f"every disease scored zero against {len(diseases):,} diseases - the "
            f"patient's terms share no informative ancestor with any of them. Either "
            f"the terms are all very general, or this hp.obo and phenotype.hpoa come "
            f"from different HPO releases.")
    top = scored[0][0]
    kept = [s for s in scored if s[0] >= top * args.min_rel_score]
    log(f"  best score {top:.3f} ({scored[0][3][:60]})")
    log(f"  diseases above {args.min_rel_score:.0%} of best: {len(kept):,}")

    loo = []
    if args.leave_one_out and len(patient) > 1:
        log("  leave-one-out (term removed -> new best disease):")
        top_id = scored[0][2]
        for drop in patient:
            sub = [q for q in patient if q != drop]
            bft = {t: max(rows_sim[q].get(t, 0.0) for q in sub) for t in all_terms}
            sc = sorted(((score_disease(sub, d["terms"], rows_sim, bft)[0], did, d["name"])
                         for did, d in diseases.items()), reverse=True)
            rank_of_top = next((i for i, x in enumerate(sc, 1) if x[1] == top_id), None)
            loo.append({"term_removed": drop, "term_name": terms[drop]["name"],
                        "new_best": sc[0][1], "new_best_name": sc[0][2],
                        "rank_of_full_best": rank_of_top})
            moved = "" if sc[0][1] == top_id else f"   [{top_id} falls to #{rank_of_top}]"
            log(f"      -{drop} {terms[drop]['name'][:26]:<26} -> "
                f"{sc[0][2][:42]}{moved}")
        if any(x["new_best"] != top_id for x in loo):
            log("    the top disease is NOT stable to dropping a single term - read the "
                "panel as one hypothesis among several, not as a ranking")

    d2g = load_disease_genes(g2p)
    coords = load_coords(args.gene_coords)
    checklist = load_checklist(args.mva_genes, args.extra_genes)

    # The ranking reads phenotype.hpoa and the genes come from genes_to_phenotype.txt, so
    # a disease can top the ranking and contribute nothing. Orphanet entries are where
    # this bites, and ORPHA:1052 is load-bearing here: BUB1 and BUB3 have no OMIM MVA
    # entry at all, so the grouped Orphanet entry is their only route into the panel.
    barren = [(f, d, n) for f, _, d, n in kept if not d2g.get(d)]
    if barren:
        log(f"  ! {len(barren)} of {len(kept)} kept diseases have no gene in "
            f"{g2p.name} and contribute nothing:")
        for f, d, n in barren[:5]:
            log(f"      {f:.3f}  {d:<14} {n[:56]}")
        if barren[0][0] == kept[0][0]:
            log("    the TOP-scoring disease is one of them - the panel is being built "
                "from weaker matches than the ranking suggests")

    best = {}
    for fwd, sym, did, name in kept:
        for gene in d2g.get(did, ()):
            cur = best.get(gene)
            if cur is None or fwd > cur["score"]:
                best[gene] = {"score": fwd, "sym": sym, "disease": did,
                              "disease_name": name, "n": cur["n"] if cur else 0}
            best[gene]["n"] += 1
    for gene in checklist:
        best.setdefault(gene, {"score": 0.0, "sym": 0.0, "disease": "",
                               "disease_name": "", "n": 0})
    log(f"  genes: {len(best):,} ({len(checklist & set(best))} of "
        f"{len(checklist)} checklist genes present)")

    log("[4/4] writing")
    rows = sorted(best.items(), key=lambda kv: (-kv[1]["score"], kv[0]))
    out = ["\t".join(["rank", "gene", "score", "score_sym", "n_diseases",
                      "best_disease", "best_disease_name", "on_checklist",
                      "chrom", "start", "end", "strand"])]
    for i, (gene, v) in enumerate(rows, 1):
        c = coords.get(gene, ("", "", "", ""))
        out.append("\t".join([str(i), gene, f"{v['score']:.4f}", f"{v['sym']:.4f}",
                              str(v["n"]), v["disease"], v["disease_name"],
                              "yes" if gene in checklist else "no", *c]))
    panel = args.outdir / "panel.tsv"
    write_atomic(panel, "\n".join(out) + "\n")

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    versions = {"hp.obo": hpo_version, "phenotype.hpoa": hpoa_version}
    params = {"min_rel_score": args.min_rel_score,
              "ranking": "Resnik MICA, best-match average, patient->disease direction",
              "ic_corpus": "diseases in phenotype.hpoa, aspect P, NOT rows dropped"}

    write_atomic(args.outdir / "panel.meta.json", json.dumps({
        "built_utc": stamp,
        "hpo_versions": versions,
        "parameters": params,
        "patient_terms": [{"id": t, "name": terms[t]["name"]} for t in patient],
        "term_notes": notes,
        "terms_not_usable": unmapped,
        "label_mismatches": mismatched,
        "n_diseases_scored": len(scored),
        "n_diseases_kept": len(kept),
        "leave_one_out": loo,
        "kept_diseases_without_genes": [{"score": round(f, 4), "id": d, "name": n}
                                        for f, d, n in barren],
        "n_genes": len(best),
        "checklist_missing_from_panel": sorted(checklist - set(best)),
        "caveats": [
            "The ranking has no null model and no p-value: the scores order diseases, "
            "they do not test them.",
            "Annotation bias is unavoidable - a well-studied syndrome carries more HPO "
            "terms and outscores an equally fitting rare one.",
            "A gene absent from this panel is NOT excluded. The panel orders the "
            "candidates coming out of the genome-wide cascade; it never filters them.",
            "Family-history terms are listed but contribute nothing to the scores. HPO "
            "disease matching treats every term as a feature of the PROBAND, so a "
            "relative's feature would pull the ranking towards diseases in which the child "
            "has that feature. They belong instead to the inheritance-model rationale, to "
            "the choice of confirmatory test, and to the differential for the relative's "
            "own phenotype - none of which this script decides.",
        ],
    }, indent=2))

    log(f"  {panel}")
    log(f"  {args.outdir / 'panel.meta.json'}  (holds the phenotype - patient material)")
    if miss := sorted(checklist - set(best)):
        log(f"  ! checklist genes with no HPO disease link: {', '.join(miss)}")
    log("\nTop 15:")
    for i, (gene, v) in enumerate(rows[:15], 1):
        flag = "*" if gene in checklist else " "
        log(f"  {i:3} {flag} {gene:<12} {v['score']:.3f}  {v['disease_name'][:52]}")


if __name__ == "__main__":
    main()
