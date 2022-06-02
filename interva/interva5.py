# -*- coding: utf-8 -*-

"""
interva.interva5
-------------------

This module contains the class for the InterVA5 algorithm.
"""

from pandas import read_csv


class InterVA5:
    """InterVA5 algorithm for assigning cause of death.

    :param va_input: Verbal Autopsy data
    :type va_input: pandas data.frame or path to CSV file
    :param hiv: likelihood of HIV as a cause of death.  Possible values are
     "H" for high (~ 1:100 deaths), "L" for low (~ 1:1000), or "V" for very
     low (~ 1:10000)
    :type hiv: string
    :param malaria: likelihood of malaria as a cause of death.  Possible values are
     "H" for high (~ 1:100 deaths), "L" for low (~ 1:1000), or "V" for very
     low (~ 1:10000)
    :type malaria: string
    """

    def __init__(self, va_input, hiv, malaria):

        self.va_input = va_input
        self.hiv = hiv
        self.malaria = malaria

    def _check_data(self):
        """Run data check."""
        pass

    def run(self):
        """Assign causes of death."""
        pass

    def get_hiv(self):
        """Get HIV parameter."""
        print(f"HIV parameter is {self.hiv}")
        return self.hiv

    def get_malaria(self):
        """Get malaria parameter."""
        print(f"Malaria parameter is {self.malaria}")
        return self.malaria

    def set_hiv(self, hiv_level):
        """Set HIV parameter."""
        # TODO: check for valid input and print message if invalid
        return self.hiv
        print(f"HIV parameter is {self.hiv}")

    def set_malaria(self, malaria_level):
        """Set malaria parameter."""
        # TODO: check for valid input and print message if invalid
        return self.malaria
        print(f"Malaria parameter is {self.malaria}")

    def get_ids(self):
        """Return pandas series of ID column in data."""
        pass

    def plot_csmf(self, top=10, file=None):
        """Plot cause-specific mortality fraction (CSMF)."""
        pass

    def get_csmf(self, top=10):
        """Print top causes in cause-specific mortality fraction (CSMF)."""
        pass

    def write_csmf(self):
        """Write cause-specific mortality fraction (CSMF) to CSV file."""
        pass

    def get_top_cod(self, top=3, include_propensities=False):
        """Get top causes of death for each individual."""
        pass

    def write_top_cod(self, top=3, include_propensities=False):
        """Write top causes of death for each individual to CSV file."""

    def get_indiv_prob(self, include_propensities=False):
        """Get individual causes of death distribution."""
        pass

    def write_indiv_prob(self, include_propensities=False):
        """Write individual cause of death distribution to CSV file."""
        pass
