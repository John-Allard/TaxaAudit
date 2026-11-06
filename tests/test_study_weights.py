import numpy as np
import pandas as pd
import pytest
from taxaaudit.model import fit_group, bootstrap_p

def frame():
    return pd.DataFrame(dict(taxid=[1]*4, name=['group']*4, rank=['family']*4,
        observed_failure=[1,1,0,0], background_apparent_failure_probability=[.1]*4,
        cross_study_detection_probability=[.8]*4, mean_sampled_members=[5]*4))

def test_constant_weights_preserve_all_results():
    f=frame(); g=np.linspace(0,1,101)
    a,s=fit_group(f,g)
    f['study_weight']=7.
    b,t=fit_group(f,g)
    assert a==b and s==t
    p=bootstrap_p(f,g,t,300,42)
    assert p==bootstrap_p(f.drop(columns='study_weight'),g,s,300,42)

def test_relative_weights_change_fitted_recurrence():
    f=frame(); g=np.linspace(0,1,101)
    baseline=fit_group(f,g)[0]['fitted_reproducible_nonmonophyly_tendency']
    f['study_weight']=[4,4,1,1]
    assert fit_group(f,g)[0]['fitted_reproducible_nonmonophyly_tendency'] > baseline

def test_weighted_bootstrap_matches_exact_null_tail():
    from itertools import product
    f=frame(); f['study_weight']=[4,2,1,1]
    grid=np.linspace(0,1,101)
    observed=fit_group(f,grid)[1]
    expected=0.
    for outcome in product([0,1], repeat=4):
        alternative=f.copy(); alternative['observed_failure']=outcome
        stat=fit_group(alternative,grid)[1]
        if stat >= observed-1e-12:
            expected += np.prod([.1 if x else .9 for x in outcome])
    assert abs(bootstrap_p(f,grid,observed,20000,42)-expected) < .01

def test_explicit_higher_rank_selection_does_not_change_default():
    from taxaaudit.taxonomy import Taxonomy
    tax=Taxonomy(parent={1:1,2:1,3:2,4:3,5:4},
        rank={1:'no rank',2:'phylum',3:'order',4:'family',5:'species'},
        name={1:'root',2:'P',3:'O',4:'F',5:'S'}, aliases={},merged={},deleted=set())
    assert tax.testable_ancestors(5)==(4,)
    tax.target_ranks=frozenset({'order','phylum'})
    assert set(tax.testable_ancestors(5))=={2,3}

@pytest.mark.parametrize('bad',[0,-1,np.nan,np.inf])
def test_invalid_weights_rejected(bad):
    f=frame(); f['study_weight']=[bad,1,1,1]
    with pytest.raises(ValueError): fit_group(f,np.linspace(0,1,11))
