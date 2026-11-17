"""Whole-read comparisons. Sites on a read are never independent samples."""
from __future__ import annotations

import collections
import math

import numpy as np


def exact_test(a, b, state_limit=250_000):
    """Exact two-sided difference of means for integer high-call counts.

    Every count must use the same positional denominator. Multiset weights
    count assignments of individual reads, including identical profiles.
    """
    if not a or not b:
        raise ValueError("Both groups must contain reads.")
    if any(int(v) != v or v < 0 for v in a+b):
        raise ValueError("Exact testing requires nonnegative integer counts.")
    values = list(map(int, a+b))
    N, n, total = len(values), len(a), sum(values)
    observed = abs(N*sum(a)-n*total)
    denominator = math.comb(N,n)
    if observed == 0:
        return {"p_value":1., "test":"exact_whole_read", "assignments":denominator}
    n = min(n, N-n)
    dp = {(0,0):1}
    for value, frequency in sorted(collections.Counter(values).items()):
        nxt = collections.defaultdict(int)
        for (k,s), weight in dp.items():
            for j in range(min(frequency,n-k)+1):
                nxt[k+j,s+j*value] += weight*math.comb(frequency,j)
                if len(nxt) > state_limit:
                    return None
        dp = nxt
    assignments = sum(weight for (k,_),weight in dp.items() if k==n)
    if assignments != denominator:
        raise RuntimeError("Exact permutation assignment count is inconsistent.")
    extreme = sum(weight for (k,s),weight in dp.items()
                  if k==n and abs(N*s-n*total) >= observed)
    return {"p_value":extreme/denominator, "test":"exact_whole_read", "assignments":denominator}


def permutation_test(a, b, seed=20261016, iterations=19_999, state_limit=250_000):
    value = exact_test(a,b,state_limit)
    if value is not None:
        return value
    rng = np.random.default_rng(seed)
    values = np.asarray(a+b,dtype=np.int64)
    n, N, total = len(a),len(values),int(values.sum())
    observed = abs(N*sum(a)-n*total)
    extreme = 0
    for _ in range(iterations):
        s = int(values[rng.choice(N,n,replace=False)].sum())
        extreme += abs(N*s-n*total) >= observed
    return {"p_value":(extreme+1)/(iterations+1), "test":"sampled_whole_read",
            "assignments":iterations, "random_seed":seed}


def bootstrap_interval(a,b,denominator,iterations,seed):
    if iterations == 0:
        return {"ci_low_pp":None,"ci_high_pp":None,"bootstrap_replicates":0}
    rng = np.random.default_rng(seed)
    a,b = np.asarray(a,dtype=float)/denominator,np.asarray(b,dtype=float)/denominator
    chunks = []
    for start in range(0,iterations,256):
        size = min(256,iterations-start)
        chunks.append(100*(rng.choice(b,(size,len(b))).mean(axis=1)-rng.choice(a,(size,len(a))).mean(axis=1)))
    low,high = np.quantile(np.concatenate(chunks),[.025,.975])
    return {"ci_low_pp":float(low),"ci_high_pp":float(high),"bootstrap_replicates":iterations}
