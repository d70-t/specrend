import specrend
import numpy as np
import numpy.testing as npt


def test_color_table():
    wvlns = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900])
    cs = specrend.CieColorTable(wvlns)
    assert cs.wvlns_idx == slice(3, 7)
    npt.assert_equal(cs.wvlns, np.array([400, 500, 600, 700]))


def test_dynamic_range():
    im = np.array([1, 2, 0])
    im = specrend.array_max_dynamic_range(im)
    npt.assert_equal(im, [0.5, 1., 0.])
