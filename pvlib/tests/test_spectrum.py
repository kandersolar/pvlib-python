import pytest
from numpy.testing import assert_allclose
import pandas as pd
import numpy as np
from pvlib import spectrum
from .conftest import DATA_DIR

SPECTRL2_TEST_DATA = DATA_DIR / 'spectrl2_example_spectra.csv'
