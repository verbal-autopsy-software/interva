# -*- coding: utf-8 -*-

import pytest
from pandas import read_csv
import pkgutil
from io import BytesIO

from interva.interva5 import InterVA5
va_data_csv = pkgutil.get_data("interva", "data/randomva5.csv")
va_data = read_csv(BytesIO(va_data_csv))


def test_get_hiv():
    iv5out = InterVA5(va_data, hiv="h", malaria="l")
    assert iv5out.get_hiv() == "h"
