import numpy as np

from badtaxafinder.model import bh_adjust


def test_bh_is_monotone_in_ranked_p_values():
    p = np.array([0.01, 0.04, 0.03, 0.8])
    q = bh_adjust(p)
    assert np.allclose(q, [0.04, 0.05333333333333334, 0.05333333333333334, 0.8])

