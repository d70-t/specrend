import specrend
import numpy as np
import numpy.testing as npt


def test_dynamic_range():
    im = np.array([1, 2, 0])
    im = specrend.array_max_dynamic_range(im)
    npt.assert_equal(im, [0.5, 1., 0.])
