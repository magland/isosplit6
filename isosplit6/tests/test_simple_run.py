from isosplit6 import isosplit6
import numpy as np
from numpy.typing import NDArray

def test_isosplit6_runs():
    """
    A simple test to check that isosplit6 runs.
    Useful for testing wheels are successfully built.
    """
    data = np.random.random((100, 100))

    result = isosplit6(data)

    assert isinstance(result, np.ndarray)