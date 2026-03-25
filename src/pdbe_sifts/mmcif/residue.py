#!/usr/bin/env python3


class Residue:
    """Represents a single residue in a PDB polypeptide chain."""

    n: int
    auth_n: str | None
    auth_ins: str | None
    oneL: str
    threeL: str
    observed: bool
    mh: int
    rtype: str | None
    oneL_original: str | None
    threeL_original: str | None
    is_chromophore: bool

    def __init__(
        self,
        n: int,
        auth_n: str | None,
        auth_ins: str | None,
        oneL: str,
        threeL: str,
        rtype: str | None,
        observed: bool,
        oneL_original: str | None,
        threeL_original: str | None,
        mh: int,
    ) -> None:
        """Initialise a Residue with all positional and type metadata.

        Args:
            n: 1-based sequence position in the mmCIF poly_seq scheme.
            auth_n: Author residue number (PDB numbering); ``None`` when
                unobserved.
            auth_ins: Insertion code; ``None`` when absent or unobserved.
            oneL: One-letter amino acid code (may be multi-character for
                chromophores).
            threeL: Three-letter residue code (e.g. ``"ALA"``).
            rtype: Annotation type from ``_struct_ref_seq_dif`` (e.g.
                ``"Engineered mutation"``, ``"Insertion"``); ``None`` for
                standard residues.
            observed: ``True`` if the residue is present in the electron
                density (observed in the structure).
            oneL_original: Original one-letter code before an engineered
                mutation or conflict; ``None`` otherwise.
            threeL_original: Original three-letter code; ``None`` otherwise.
            mh: Microheterogeneity index (1 for the primary conformer, >1
                for alternate conformers at the same position).
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

    def __repr__(self) -> str:
        """Return a compact string representation of the residue.

        Returns:
            String of the form
            ``"<n> (<auth_n><auth_ins>) [<oneL>|<threeL>] T: <rtype> O: <observed> MH: <mh>"``
            useful for logging and debugging.
        """
        return f"{self.n} ({self.auth_n}{self.auth_ins}) [{self.oneL}|{self.threeL}] T: {self.rtype} O: {self.observed} MH: {self.mh}"
