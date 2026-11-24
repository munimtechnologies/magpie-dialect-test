# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""magpie-dialect-test: pre-registered re-test of group-specific call sequences in
Western Australian magpie fledglings (Mason, King and Ridley 2026, Proc. R. Soc. B 293:20251620).

Tests (PLAN.md v0 as amended by AMENDMENT-1.md, 2026-11-24):
  T0  the authors' binomial GLMMs fl01 ~ group01 / exgroup01 / 1 + (1|ID), adaptive
      Gauss-Hermite maximum likelihood; published-style decision.
  T1  own-group overlap excess D, null = fledglings relabeled between groups (counts kept).
  T2  order beyond content among multi-order bags, same null.
  T3  donor worlds: each fledgling gets a different group's real records; how often T0's
      decision and T1 say "own group".
  T4  social contacts against effort-adjusted repertoire size, Spearman, permutation null.
  Power for T1 and T2 with planted effects.
Run: python magpie_test.py --data-dir DIR --out-dir OUT
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

SEED = 20261124
ALPHA = 0.05
GH_NODES = 25
RELABEL_WORLDS = 9999
DONOR_WORLDS = 200
DONOR_RELABELS = 999
POWER_WORLDS = 50
POWER_RELABELS = 999
T4_PERMS = 99999
MISCAL_LIMIT = 0.10
PRIMARY = "COMBINATION_repertoire.csv"
SECONDARY = "COMBINATION_repertoire_single_occurrences_removed.csv"
FREQ = "Frequency_diversity_data.csv"


# ---------------------------------------------------------------- loading

def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def norm_type(s):
    return " ".join(s.split())


def build(rows):
    """Fledgling x type presence matrices from the authors' long-format rows."""
    if not rows:
        raise ValueError("no rows")
    for k in ("ID", "group", "call", "fl01", "group01", "exgroup01"):
        if k not in rows[0]:
            raise ValueError("missing column %s" % k)
    ids = sorted({r["ID"].strip() for r in rows})
    types = sorted({norm_type(r["call"]) for r in rows})
    fi = {f: i for i, f in enumerate(ids)}
    ti = {t: j for j, t in enumerate(types)}
    n, m = len(ids), len(types)
    F = np.zeros((n, m), dtype=np.int64)
    G = np.zeros((n, m), dtype=np.int64)
    E = np.zeros((n, m), dtype=np.int64)
    seen = np.zeros((n, m), dtype=bool)
    group_of = {}
    for r in rows:
        f, g, t = r["ID"].strip(), r["group"].strip(), norm_type(r["call"])
        if group_of.setdefault(f, g) != g:
            raise ValueError("fledgling %s listed in two groups" % f)
        vals = []
        for k in ("fl01", "group01", "exgroup01"):
            v = r[k].strip()
            if v not in ("0", "1"):
                raise ValueError("value %r in %s for %s / %s" % (v, k, f, t))
            vals.append(int(v))
        i, j = fi[f], ti[t]
        if seen[i, j]:
            raise ValueError("duplicate row %s / %s" % (f, t))
        seen[i, j] = True
        F[i, j], G[i, j], E[i, j] = vals
    groups = sorted(set(group_of.values()))
    gi = {g: k for k, g in enumerate(groups)}
    gidx = np.array([gi[group_of[f]] for f in ids], dtype=np.int64)
    return {"ids": ids, "types": types, "groups": groups, "gidx": gidx,
            "F": F, "G": G, "E": E, "seen": seen, "rows": len(rows)}


def group_records(G, gidx, n_groups):
    """R[g] = union of the follow records of group g's fledglings."""
    R = np.zeros((n_groups, G.shape[1]), dtype=np.int64)
    for f in range(G.shape[0]):
        R[gidx[f]] = np.maximum(R[gidx[f]], G[f])
    return R


def loader_checks(d):
    F, G, E, gidx = d["F"], d["G"], d["E"], d["gidx"]
    out = {"rows": d["rows"], "fledglings": len(d["ids"]), "types": len(d["types"]),
           "groups": len(d["groups"]), "all_cells_present": bool(d["seen"].all()),
           "fledglings_per_group": {g: int(np.sum(gidx == k)) for k, g in enumerate(d["groups"])}}
    same = {}
    for k, g in enumerate(d["groups"]):
        mem = np.flatnonzero(gidx == k)
        if len(mem) > 1:
            same[g] = bool(all(np.array_equal(G[mem[0]], G[h]) for h in mem[1:]))
    out["multi_fledgling_group_records_identical"] = same
    agree = []
    for f in range(F.shape[0]):
        other = G[gidx != gidx[f]].max(axis=0)
        agree.append(float(np.mean(other == E[f])))
    out["exgroup_equals_union_of_other_groups_share"] = agree
    inside = []
    for f in range(F.shape[0]):
        nf = F[f].sum()
        inside.append(float((F[f] * G[f]).sum() / nf) if nf else None)
    out["share_of_F_in_own_follow_record"] = inside
    out["repertoire_sizes"] = [int(x) for x in F.sum(axis=1)]
    out["follow_record_sizes"] = [int(x) for x in G.sum(axis=1)]
    return out


# ---------------------------------------------------------------- statistics

def perm_p(obs, null):
    null = np.asarray([x for x in null if x is not None and np.isfinite(x)], dtype=float)
    if obs is None or not np.isfinite(obs) or null.size == 0:
        return None
    return float((1 + np.sum(null >= obs - 1e-12)) / (1 + null.size))


def overlap_matrix(F, R):
    """M[f, g] = |F_f and R_g| / |F_f| (rows with empty F_f must be removed first)."""
    sizes = F.sum(axis=1).astype(float)
    if np.any(sizes == 0):
        raise ValueError("empty repertoire")
    return (F @ R.T) / sizes[:, None]


def d_stat(M, assign):
    """Mean over fledglings of own-group overlap minus mean overlap with the other groups."""
    n, k = M.shape
    own = M[np.arange(n), assign]
    other = (M.sum(axis=1) - own) / (k - 1)
    return float(np.mean(own - other))


def relabel_null(stat_fn, assign, worlds, rng):
    return [stat_fn(rng.permutation(assign)) for _ in range(worlds)]


def t1_test(F, R, assign, worlds, rng):
    M = overlap_matrix(F, R)
    obs = d_stat(M, assign)
    null = relabel_null(lambda a: d_stat(M, a), assign, worlds, rng)
    p = perm_p(obs, null)
    return {"D": obs, "worlds": worlds, "p": p, "null_mean": float(np.mean(null)),
            "null_q95": float(np.quantile(null, 0.95)),
            "reading": "OWN_GROUP_EXCESS" if (obs > 0 and p is not None and p <= ALPHA)
            else "NO_OWN_GROUP_EXCESS",
            "own_overlap": [float(x) for x in M[np.arange(len(assign)), assign]],
            "mean_other_overlap": [float(x) for x in
                                   (M.sum(axis=1) - M[np.arange(len(assign)), assign])
                                   / (M.shape[1] - 1)]}


def bag(t):
    return tuple(sorted(t.split()))


def bag_info(types):
    bags = sorted({bag(t) for t in types})
    bi = {b: k for k, b in enumerate(bags)}
    bag_id = np.array([bi[bag(t)] for t in types], dtype=np.int64)
    orders = defaultdict(set)
    for t in types:
        orders[bag(t)].add(t)
    multi = np.array([len(orders[bag(t)]) >= 2 for t in types], dtype=bool)
    return bag_id, multi


def order_matrices(F, R, bag_id, multi):
    """A[f,g]: f's multi-order types present in R_g in f's order; B[f,g]: same bag present."""
    n_groups = R.shape[0]
    Bpres = np.zeros_like(R)
    for g in range(n_groups):
        present = set(bag_id[R[g] == 1].tolist())
        Bpres[g] = np.array([b in present for b in bag_id], dtype=np.int64)
    Fm = F * multi.astype(np.int64)[None, :]
    return Fm @ R.T, Fm @ Bpres.T


def d_ord(A, B, assign):
    n, k = A.shape
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = np.where(B > 0, A / np.where(B > 0, B, 1), np.nan)
    terms = []
    for f in range(n):
        own = Q[f, assign[f]]
        others = np.array([Q[f, g] for g in range(k) if g != assign[f]])
        others = others[np.isfinite(others)]
        if not np.isfinite(own) or others.size == 0:
            continue
        terms.append(own - others.mean())
    if not terms:
        return None, 0
    return float(np.mean(terms)), len(terms)


def t2_test(F, R, assign, types, worlds, rng):
    bag_id, multi = bag_info(types)
    produced_multi = int(np.sum(multi & (F.sum(axis=0) > 0)))
    A, B = order_matrices(F, R, bag_id, multi)
    obs, entered = d_ord(A, B, assign)
    null = [d_ord(A, B, rng.permutation(assign))[0] for _ in range(worlds)]
    p = perm_p(obs, null)
    if produced_multi < 10 or entered < 5:
        reading = "UNINFORMATIVE"
    elif obs is not None and obs > 0 and p is not None and p <= ALPHA:
        reading = "ORDER_CARRIES_GROUP_SIGNAL"
    else:
        reading = "ORDER_ADDS_NOTHING_MEASURABLE"
    vals = [x for x in null if x is not None]
    return {"D_ord": obs, "fledglings_entered": entered,
            "multi_order_types_in_universe": int(multi.sum()),
            "multi_order_types_produced": produced_multi, "worlds": len(vals), "p": p,
            "null_mean": float(np.mean(vals)) if vals else None, "reading": reading}


# ---------------------------------------------------------------- GLMM (T0)

def _expit(x):
    return 0.5 * (1.0 + np.tanh(0.5 * x))


def glmm_nll(theta, y, x, onehot, cluster, nodes, logw):
    """Negative marginal log-likelihood of a binomial random-intercept model, adaptive GH."""
    b0 = theta[0]
    b1 = theta[1] if x is not None else 0.0
    s = math.exp(theta[-1])
    eta = b0 + (b1 * x if x is not None else 0.0)
    C = onehot.shape[0]
    u = np.zeros(C)
    for _ in range(30):
        p = _expit(eta + u[cluster])
        grad = onehot @ (y - p) - u / s ** 2
        hess = -(onehot @ (p * (1.0 - p))) - 1.0 / s ** 2
        step = np.clip(grad / hess, -5.0, 5.0)
        u = u - step
        if np.max(np.abs(step)) < 1e-12:
            break
    p = _expit(eta + u[cluster])
    hess = -(onehot @ (p * (1.0 - p))) - 1.0 / s ** 2
    tau = 1.0 / np.sqrt(-hess)
    U = u[:, None] + math.sqrt(2.0) * tau[:, None] * nodes[None, :]       # (C, K)
    lin = eta[:, None] + U[cluster] if np.ndim(eta) else eta + U[cluster]
    ll = y[:, None] * lin - np.logaddexp(0.0, lin)                        # (n, K)
    S = onehot @ ll - U ** 2 / (2.0 * s ** 2)                             # (C, K)
    T = S + logw[None, :] + nodes[None, :] ** 2
    mx = T.max(axis=1)
    lse = mx + np.log(np.exp(T - mx[:, None]).sum(axis=1))
    logint = -0.5 * math.log(2.0 * math.pi * s ** 2) + 0.5 * math.log(2.0) + np.log(tau) + lse
    return float(-logint.sum())


def num_hessian(fn, theta, h=1e-4):
    k = len(theta)
    H = np.zeros((k, k))
    for i in range(k):
        for j in range(i, k):
            def f(di, dj):
                t = np.array(theta, dtype=float)
                t[i] += di
                t[j] += dj
                return fn(t)
            H[i, j] = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
            H[j, i] = H[i, j]
    return H


def glmm_fit(y, x, cluster, n_clusters):
    from scipy.optimize import minimize
    y = np.asarray(y, dtype=float)
    x = None if x is None else np.asarray(x, dtype=float)
    cluster = np.asarray(cluster, dtype=np.int64)
    onehot = np.zeros((n_clusters, len(y)))
    onehot[cluster, np.arange(len(y))] = 1.0
    z, w = np.polynomial.hermite.hermgauss(GH_NODES)
    logw = np.log(w)

    def fn(t):
        return glmm_nll(t, y, x, onehot, cluster, z, logw)

    ybar = min(max(y.mean(), 1e-3), 1 - 1e-3)
    x0 = [math.log(ybar / (1 - ybar))] + ([0.0] if x is not None else []) + [math.log(0.5)]
    bounds = [(-20, 20)] + ([(-20, 20)] if x is not None else []) + [(-6, 3)]
    res = minimize(fn, x0, method="L-BFGS-B", bounds=bounds)
    nll = float(res.fun)
    k = len(x0)
    n = len(y)
    aicc = 2 * k + 2 * nll + 2 * k * (k + 1) / (n - k - 1)
    out = {"params": [float(v) for v in res.x], "nll": nll, "aicc": aicc,
           "converged": bool(res.success), "b1": None, "se": None, "z": None, "p": None}
    if x is not None:
        out["b1"] = float(res.x[1])
        try:
            cov = np.linalg.inv(num_hessian(fn, res.x))
            var = cov[1, 1]
            if np.isfinite(var) and var > 0:
                se = math.sqrt(var)
                zz = res.x[1] / se
                out.update({"se": se, "z": float(zz), "p": float(math.erfc(abs(zz) / math.sqrt(2)))})
        except np.linalg.LinAlgError:
            pass
    return out


def t0_decision(F, Gx, Ex, seen):
    """The authors' three models on the fledgling x type cells present in the file."""
    ii, jj = np.nonzero(seen)
    y = F[ii, jj]
    base = glmm_fit(y, None, ii, F.shape[0])
    pop = glmm_fit(y, Ex[ii, jj], ii, F.shape[0])
    grp = glmm_fit(y, Gx[ii, jj], ii, F.shape[0])
    fires = bool(grp["b1"] is not None and grp["b1"] > 0 and grp["p"] is not None
                 and grp["p"] <= ALPHA and grp["aicc"] < pop["aicc"] and grp["aicc"] < base["aicc"])
    return {"base": base, "pop": pop, "group": grp, "decision_fires": fires,
            "cells": int(len(y))}


# ---------------------------------------------------------------- T3 donor worlds

def donor_map(gidx, rng, tries=200000):
    n = len(gidx)
    for _ in range(tries):
        p = rng.permutation(n)
        if np.all(gidx[p] != gidx):
            return p
    raise RuntimeError("no donor map found")


def t3_calibration(d, R, worlds, relabels, rng, deadline, run_glmm=True):
    F, G, E, gidx, seen = d["F"], d["G"], d["E"], d["gidx"], d["seen"]
    keep = F.sum(axis=1) > 0
    t0_hits = t1_hits = done = 0
    for _ in range(worlds):
        if time.time() > deadline:
            break
        p = donor_map(gidx, rng)
        if run_glmm:
            r0 = t0_decision(F, G[p], E[p], seen)
            t0_hits += int(r0["decision_fires"])
        assign = gidx[p][keep]
        r1 = t1_test(F[keep], R, assign, relabels, rng)
        t1_hits += int(r1["reading"] == "OWN_GROUP_EXCESS")
        done += 1
    out = {"worlds": done, "t1_hits": t1_hits,
           "t1_rate": (t1_hits / done) if done else None}
    if run_glmm:
        out.update({"t0_hits": t0_hits, "t0_rate": (t0_hits / done) if done else None})
        out["t0_reading"] = ("MISCALIBRATED" if done and t0_hits / done > MISCAL_LIMIT
                             else "CALIBRATED") if done else "INCOMPLETE"
    return out


# ---------------------------------------------------------------- power

def plant_t1(F, R, assign, q, rng):
    Fp = F.copy()
    for f in range(F.shape[0]):
        own = R[assign[f]]
        for j in np.flatnonzero((F[f] == 1) & (own == 0)):
            if rng.random() < q:
                cand = np.flatnonzero((own == 1) & (Fp[f] == 0))
                if cand.size:
                    Fp[f, j] = 0
                    Fp[f, rng.choice(cand)] = 1
    return Fp


def plant_t2(F, R, assign, bag_id, multi, rng):
    Fp = F.copy()
    for f in range(F.shape[0]):
        own = R[assign[f]]
        for j in np.flatnonzero((F[f] == 1) & multi & (own == 0)):
            cand = np.flatnonzero((bag_id == bag_id[j]) & (own == 1) & (Fp[f] == 0))
            if cand.size and rng.random() < 0.5:
                Fp[f, j] = 0
                Fp[f, rng.choice(cand)] = 1
    return Fp


def power(F, R, assign, types, q, rng, deadline):
    bag_id, multi = bag_info(types)
    h1 = h2 = d1 = d2 = 0
    for _ in range(POWER_WORLDS):
        if time.time() > deadline:
            break
        r = t1_test(plant_t1(F, R, assign, q, rng), R, assign, POWER_RELABELS, rng)
        h1 += int(r["reading"] == "OWN_GROUP_EXCESS")
        d1 += 1
        r2 = t2_test(plant_t2(F, R, assign, bag_id, multi, rng), R, assign, types,
                     POWER_RELABELS, rng)
        h2 += int(r2["reading"] == "ORDER_CARRIES_GROUP_SIGNAL")
        d2 += 1
    return {"q_t1": q, "t1": {"worlds": d1, "hits": h1, "power": (h1 / d1) if d1 else None},
            "t2": {"worlds": d2, "hits": h2, "power": (h2 / d2) if d2 else None}}


# ---------------------------------------------------------------- T4

def rankdata(a):
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a))
    sa = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = math.sqrt(float(ra @ ra) * float(rb @ rb))
    return float(ra @ rb / den) if den > 0 else None


def spearman_perm_p(contacts, target, perms, rng):
    obs = spearman(contacts, target)
    if obs is None:
        return None, None
    rc = rankdata(contacts)
    rt = rankdata(target)
    rc, rt = rc - rc.mean(), rt - rt.mean()
    den = math.sqrt(float(rc @ rc) * float(rt @ rt))
    hits, done = 0, 0
    while done < perms:
        b = min(10000, perms - done)
        idx = rng.permuted(np.tile(np.arange(len(rc)), (b, 1)), axis=1)
        vals = (rc[idx] @ rt) / den
        hits += int(np.sum(vals >= obs - 1e-12))
        done += b
    return obs, float((1 + hits) / (1 + done))


def t4_test(freq_rows, ids, sizes, rng):
    by = defaultdict(list)
    for r in freq_rows:
        by[r["ID"].strip()].append(r)
    missing = [f for f in ids if f not in by]
    contacts, effort, size, used, na, nonconst = [], [], [], [], 0, []
    for f, s in zip(ids, sizes):
        if f not in by:
            continue
        cs = {r["total_soc_contacts"].strip() for r in by[f]}
        if len(cs) != 1:
            nonconst.append(f)
        try:
            c = float(sorted(cs)[0])
        except ValueError:
            continue
        e = 0.0
        for r in by[f]:
            try:
                e += float(r["combo_freq"])
            except ValueError:
                na += 1
        contacts.append(c)
        effort.append(e)
        size.append(float(s))
        used.append(f)
    out = {"fledglings_used": used, "missing_in_freq_file": missing,
           "contacts_not_constant": nonconst, "combo_freq_unparsed": na,
           "contacts": contacts, "effort": effort, "size": size}
    if len(used) < 8:
        out["reading"] = "NOT_RUN"
        return out
    ls, le = np.log1p(size), np.log1p(effort)
    X = np.column_stack([np.ones(len(le)), le])
    beta, *_ = np.linalg.lstsq(X, ls, rcond=None)
    resid = ls - X @ beta
    rho, p = spearman_perm_p(contacts, resid, T4_PERMS, rng)
    out.update({"rho_adjusted": rho, "p": p, "perms": T4_PERMS,
                "rho_raw_contacts_size": spearman(contacts, size),
                "rho_effort_size": spearman(effort, size),
                "reading": "CONTACTS_EFFECT_SURVIVES_EFFORT"
                if (rho is not None and rho > 0 and p is not None and p <= ALPHA)
                else "CONTACTS_EFFECT_NOT_SEPARABLE_FROM_EFFORT"})
    return out


# ---------------------------------------------------------------- main

def analyse(d, rng, deadline, with_t2_t4=True, freq_rows=None):
    keep = d["F"].sum(axis=1) > 0
    R = group_records(d["G"], d["gidx"], len(d["groups"]))
    res = {"loader": loader_checks(d), "dropped_empty_fledglings":
           [f for f, k in zip(d["ids"], keep) if not k]}
    Fk, ak = d["F"][keep], d["gidx"][keep]
    res["T0"] = t0_decision(d["F"], d["G"], d["E"], d["seen"])
    print("T0 done", res["T0"]["decision_fires"], flush=True)
    res["T1"] = t1_test(Fk, R, ak, RELABEL_WORLDS, rng)
    print("T1 done p=%s" % res["T1"]["p"], flush=True)
    if with_t2_t4:
        res["T2"] = t2_test(Fk, R, ak, d["types"], RELABEL_WORLDS, rng)
        print("T2 done %s" % res["T2"]["reading"], flush=True)
    res["T3"] = t3_calibration(d, R, DONOR_WORLDS, DONOR_RELABELS, rng, deadline)
    print("T3 done %s" % res["T3"], flush=True)
    if with_t2_t4:
        q = max(0.10, res["T1"]["D"])
        res["power"] = power(Fk, R, ak, d["types"], q, rng, deadline)
        print("power done %s" % res["power"], flush=True)
        ids_k = [f for f, k in zip(d["ids"], keep) if k]
        res["T4"] = t4_test(freq_rows, ids_k, Fk.sum(axis=1).tolist(), rng)
        print("T4 done %s" % res["T4"].get("reading"), flush=True)
    return res


def predictions(res):
    t2 = res["T2"]["reading"]
    out = {"Q1_T1_own_group_excess": res["T1"]["reading"] == "OWN_GROUP_EXCESS",
           "Q2_T2_order_carries_group_signal": t2 == "ORDER_CARRIES_GROUP_SIGNAL",
           "Q3_T0_decision_miscalibrated": res["T3"].get("t0_reading") == "MISCALIBRATED",
           "Q4_T4_contacts_survive_effort": res["T4"].get("reading") == "CONTACTS_EFFECT_SURVIVES_EFFORT",
           "Q5_T2_uninformative": t2 == "UNINFORMATIVE",
           "Q6_T0_fires_in_over_half_of_donor_worlds":
               bool(res["T3"].get("t0_rate") is not None and res["T3"]["t0_rate"] > 0.5)}
    priors = {"Q1_T1_own_group_excess": 0.60, "Q2_T2_order_carries_group_signal": 0.25,
              "Q3_T0_decision_miscalibrated": 0.40, "Q4_T4_contacts_survive_effort": 0.35,
              "Q5_T2_uninformative": 0.30, "Q6_T0_fires_in_over_half_of_donor_worlds": 0.55}
    brier = {k: (priors[k] - (1.0 if v else 0.0)) ** 2 for k, v in out.items()}
    return {"outcomes": out, "priors": priors, "brier": brier,
            "mean_brier": float(np.mean(list(brier.values())))}


def write_md(res, path):
    p = res["primary"]
    L = ["# Magpie group-sequence re-test: results", "",
         "Made with AI agents working in MV2's name (MV2, Munim, Inc.); not yet checked by a human expert.", "",
         "## Primary file (%s)" % PRIMARY, "",
         "- T0 group model: b1 %s, z %s, p %s; AICc group %.2f, pop %.2f, base %.2f; published-style decision fires: %s"
         % (p["T0"]["group"]["b1"], p["T0"]["group"]["z"], p["T0"]["group"]["p"],
            p["T0"]["group"]["aicc"], p["T0"]["pop"]["aicc"], p["T0"]["base"]["aicc"],
            p["T0"]["decision_fires"]),
         "- T1: D %.4f, null mean %.4f, null q95 %.4f, p %s -> %s"
         % (p["T1"]["D"], p["T1"]["null_mean"], p["T1"]["null_q95"], p["T1"]["p"], p["T1"]["reading"]),
         "- T2: D_ord %s over %s fledglings, p %s, multi-order types produced %s -> %s"
         % (p["T2"]["D_ord"], p["T2"]["fledglings_entered"], p["T2"]["p"],
            p["T2"]["multi_order_types_produced"], p["T2"]["reading"]),
         "- T3 donor worlds: %s; T0 decision fired in %s of %s (%s); T1 fired in %s"
         % (p["T3"]["worlds"], p["T3"].get("t0_hits"), p["T3"]["worlds"], p["T3"].get("t0_reading"),
            p["T3"]["t1_hits"]),
         "- Power: %s" % json.dumps(p["power"]),
         "- T4: rho adjusted %s, p %s, raw rho %s, effort-size rho %s -> %s"
         % (p["T4"].get("rho_adjusted"), p["T4"].get("p"), p["T4"].get("rho_raw_contacts_size"),
            p["T4"].get("rho_effort_size"), p["T4"].get("reading")),
         "", "## Secondary file (%s)" % SECONDARY, ""]
    s = res.get("secondary")
    if s:
        L += ["- T0 decision fires: %s (b1 %s, p %s)" % (s["T0"]["decision_fires"], s["T0"]["group"]["b1"],
                                                        s["T0"]["group"]["p"]),
              "- T1: D %.4f, p %s -> %s" % (s["T1"]["D"], s["T1"]["p"], s["T1"]["reading"]),
              "- T3: T0 fired in %s of %s donor worlds; T1 in %s" % (s["T3"].get("t0_hits"),
                                                                   s["T3"]["worlds"], s["T3"]["t1_hits"])]
    else:
        L.append("- not run: %s" % res.get("secondary_error"))
    L += ["", "## Predictions", "", "```", json.dumps(res["predictions"], indent=1), "```", "",
          "## Loader checks", "", "```", json.dumps(p["loader"], indent=1), "```", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--budget-min", type=float, default=150.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    t_start = time.time()
    deadline = t_start + 60 * a.budget_min
    rng = np.random.default_rng(SEED)
    d = build(read_csv(os.path.join(a.data_dir, PRIMARY)))
    freq = read_csv(os.path.join(a.data_dir, FREQ))
    res = {"paper": "doi:10.1098/rspb.2025.1620", "data": "doi:10.6084/m9.figshare.28251398.v2",
           "seed": SEED, "settings": {"relabel_worlds": RELABEL_WORLDS, "donor_worlds": DONOR_WORLDS,
                                      "donor_relabels": DONOR_RELABELS, "power_worlds": POWER_WORLDS,
                                      "t4_perms": T4_PERMS, "gh_nodes": GH_NODES}}
    fr = os.path.join(a.data_dir, "fetch_report.json")
    if os.path.exists(fr):
        with open(fr, encoding="utf-8") as fh:
            res["fetch_report"] = json.load(fh)
    res["primary"] = analyse(d, rng, deadline, True, freq)
    try:
        d2 = build(read_csv(os.path.join(a.data_dir, SECONDARY)))
        res["secondary"] = analyse(d2, rng, deadline, False)
    except Exception as exc:  # noqa: BLE001  secondary is reported, never gating
        res["secondary"] = None
        res["secondary_error"] = repr(exc)
    p = res["primary"]
    res["predictions"] = predictions(p)
    res["elapsed_s"] = time.time() - t_start
    with open(os.path.join(a.out_dir, "magpie_results.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    write_md(res, os.path.join(a.out_dir, "RESULTS-MAGPIE.md"))
    print("done in %.0fs" % res["elapsed_s"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
