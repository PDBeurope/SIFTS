#!/usr/bin/env python3


class Residue:
    """Docstring for Residue."""

    def __init__(
        self,
        n,
        auth_n,
        auth_ins,
        oneL,
        threeL,
        rtype,
        observed,
        oneL_original,
        threeL_original,
        mh,
    ):
        """tt

        @param n TODO
        @param one_code TODO

        """
        self.n = n
        self.auth_n = auth_n
        self.auth_ins = auth_ins
        self.oneL = oneL
        self.threeL = threeL
        self.observed = observed
        self.mh = mh
        self.rtype = rtype
        self.oneL_original = oneL_original
        self.threeL_original = threeL_original

        self.is_chromophore = self.rtype == "Chromophore"

    def __repr__(self):
        return f"{self.n} ({self.auth_n}{self.auth_ins}) [{self.oneL}|{self.threeL}] T: {self.rtype} O: {self.observed} MH: {self.mh}"
