from fractions import Fraction
from itertools import combinations
import math

import pytest

from metaphaser.stats import exact_test, permutation_test, bootstrap_interval


def brute(a,b):
    values=a+b
    N,n,total=len(values),len(a),sum(values)
    observed=abs(N*sum(a)-n*total)
    return sum(abs(N*sum(values[i] for i in ids)-n*total)>=observed
               for ids in combinations(range(N),n))/math.comb(N,n)


@pytest.mark.parametrize("a,b",[([0,0,0],[0,0,0]),([1,1,2],[0,0,0]),
                              ([1,1,1,0],[0,0,0]),([0,2,5],[1,1,3,6]),([1],[0,2,5])])
def test_exact_probability_matches_individual_read_enumeration(a,b):
    assert exact_test(a,b)["p_value"]==pytest.approx(brute(a,b))
    assert exact_test(a,b)["assignments"]==math.comb(len(a+b),len(a))


def test_saved_setun06_probability():
    result=exact_test([0]*10+[1]*11+[2],[0]*21)
    assert result["p_value"]==pytest.approx(6.132053778111634e-5,rel=1e-13)
    assert result["assignments"]==math.comb(43,22)


def test_ties_are_included_and_group_swap_does_not_change_p():
    a,b=[0,1,2,2],[0,2,2]
    assert exact_test(a,b)["p_value"]==exact_test(b,a)["p_value"]


def test_fallback_is_labelled_reproducible_and_nonzero():
    r=permutation_test([2,2,1],[0,0,0],seed=1,iterations=99,state_limit=1)
    assert r==permutation_test([2,2,1],[0,0,0],seed=1,iterations=99,state_limit=1)
    assert r["test"]=="sampled_whole_read"
    assert r["p_value"]>=.01


def test_bootstrap_keeps_fixed_profiles_and_units():
    r=bootstrap_interval([3,3,3],[0,0,0],3,100,42)
    assert r["ci_low_pp"]==r["ci_high_pp"]==-100
    assert r["bootstrap_replicates"]==100
