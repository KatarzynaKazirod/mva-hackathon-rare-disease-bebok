#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Minimal Open Targets GraphQL client - standard library only.

The original analysis reached Open Targets through an MCP connector. That
connector is specific to the environment it ran in, so every query it made is
re-expressed here against the public endpoint

    https://api.platform.opentargets.org/api/v4/graphql

so that anyone can re-run the analysis with nothing but a Python interpreter.

Two schema differences between the connector and the public v4 endpoint are
worth knowing before you read the pipeline:

  * `Drug` on the public endpoint exposes `id`, `name`, `tradeNames`,
    `maximumClinicalStage` and `mechanismsOfAction` - and nothing else. There
    is no `hasBeenWithdrawn` and no `blackBoxWarning` field, so those two
    columns cannot be regenerated from here and are carried as curated values.
  * The drug connection is `target.drugAndClinicalCandidates` and takes no
    arguments; its rows are `ClinicalTargetFromTarget`, which carries
    `maxClinicalStage` alongside the nested `drug`.

Open Targets releases are versioned. Counts in this repository were taken from
the release live in October 2026; a later release can legitimately return
different numbers. Always record what `meta()` returns next to any count.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request

ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"


class OpenTargetsError(RuntimeError):
    """A GraphQL error, or an HTTP failure that outlived the retries."""


def query(gql: str, variables: dict | None = None, *, retries: int = 4,
          timeout: int = 90, endpoint: str = ENDPOINT) -> dict:
    """POST a GraphQL query and return `data`, retrying transient failures.

    GraphQL reports schema and validation problems as HTTP 400 with the detail
    in the body, so the body is raised rather than swallowed - a silent empty
    result is how a census quietly becomes wrong.
    """
    payload = json.dumps({"query": gql, "variables": variables or {}}).encode()
    last = None
    for attempt in range(retries):
        req = urllib.request.Request(
            endpoint, data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                body = json.load(fh)
            if body.get("errors"):
                raise OpenTargetsError(json.dumps(body["errors"])[:1000])
            return body["data"]
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode(errors="replace")[:1000]
            if exc.code == 400:
                raise OpenTargetsError(f"HTTP 400 - query rejected: {detail}") from exc
            last = f"HTTP {exc.code}: {detail}"
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = f"{type(exc).__name__}: {exc}"
        if attempt < retries - 1:
            time.sleep(2 ** attempt)
    raise OpenTargetsError(f"giving up after {retries} attempts - {last}")


def meta(**kw) -> dict:
    """Data version and API version, to be recorded beside any count."""
    d = query("{ meta { name apiVersion{x y z} dataVersion{year month iteration} } }", **kw)
    m = d["meta"]
    api = m.get("apiVersion") or {}
    dv = m.get("dataVersion") or {}
    return {
        "api_version": ".".join(str(api[k]) for k in ("x", "y", "z") if api.get(k) is not None),
        "data_version": ".".join(str(dv[k]) for k in ("year", "month", "iteration")
                                 if dv.get(k) is not None),
        "name": m.get("name"),
    }


def resolve_symbol(symbol: str, **kw) -> str | None:
    """Map an approved symbol to an Ensembl gene id, or None if unresolved.

    Search returns near-misses as well as exact hits, so the approved symbol is
    compared explicitly. Accepting the first hit is how `MT` becomes `MCAT` and
    a junk token silently enters the node set.
    """
    gql = """query($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:10})
             { hits{ id object{ ... on Target{ approvedSymbol } } } } }"""
    for hit in query(gql, {"q": symbol}, **kw)["search"]["hits"]:
        obj = hit.get("object") or {}
        if (obj.get("approvedSymbol") or "").upper() == symbol.upper():
            return hit["id"]
    return None


TARGET_DRUGS = """query($id:String!){
  target(ensemblId:$id){
    id approvedSymbol
    drugAndClinicalCandidates{
      count
      rows{
        maxClinicalStage
        drug{
          id name maximumClinicalStage
          mechanismsOfAction{ rows{ actionType mechanismOfAction targets{ approvedSymbol } } }
        }
      }
    }
  }
}"""


def target_drugs(ensembl_id: str, **kw) -> dict | None:
    """Every drug and clinical candidate recorded against one target.

    Returns None when the id is not in the release, which is different from a
    target that is present and simply has no drugs - the census must not
    collapse those two cases into one zero.
    """
    return query(TARGET_DRUGS, {"id": ensembl_id}, **kw)["target"]


def unique_drugs(target_record: dict, symbol: str | None = None) -> dict[str, dict]:
    """Collapse drug rows to one entry per ChEMBL id.

    One drug can appear in several rows (one per indication), so counting rows
    overstates the number of compounds. Action types are kept only from
    mechanisms that actually name `symbol`: a multi-target drug otherwise
    donates the action type it has on some unrelated protein, which would flip
    the direction gate on evidence about a different gene.
    """
    out: dict[str, dict] = {}
    for row in (target_record.get("drugAndClinicalCandidates") or {}).get("rows") or []:
        drug = row.get("drug")
        if not drug:
            continue
        entry = out.setdefault(drug["id"], {
            "drug_id": drug["id"],
            "drug": drug.get("name"),
            "max_stage_drug": drug.get("maximumClinicalStage"),
            "action_types": set(),
            "moa": set(),
        })
        for mech in (drug.get("mechanismsOfAction") or {}).get("rows") or []:
            named = {t.get("approvedSymbol") for t in (mech.get("targets") or [])}
            if symbol is not None and symbol not in named:
                continue
            if mech.get("actionType"):
                entry["action_types"].add(mech["actionType"])
            if mech.get("mechanismOfAction"):
                entry["moa"].add(mech["mechanismOfAction"])
    return out


if __name__ == "__main__":
    info = meta()
    print(f"endpoint   : {ENDPOINT}")
    print(f"api        : {info['api_version']}")
    print(f"data release: {info['data_version']}")
    for sym in sys.argv[1:]:
        eid = resolve_symbol(sym)
        if eid is None:
            print(f"{sym:10s} UNRESOLVED")
            continue
        rec = target_drugs(eid)
        drugs = unique_drugs(rec, rec["approvedSymbol"])
        appr = sum(1 for d in drugs.values() if d["max_stage_drug"] == "APPROVAL")
        print(f"{sym:10s} {eid}  drugs={len(drugs):3d}  approved={appr}")
