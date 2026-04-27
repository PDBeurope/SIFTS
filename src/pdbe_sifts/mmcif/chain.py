#!/usr/bin/env python3
from itertools import groupby
from multiprocessing.dummy import Pool
from operator import itemgetter

import tqdm
from Bio.Seq import Seq

from pdbe_sifts.base.utils import get_cpu_count
from pdbe_sifts.mmcif import mmcif_helper
from pdbe_sifts.segments_generation.connectivity.process_connectivity import (
    ConnectivityCheck,
)

from ..taxonomy_fix_pkl import TaxonomyFix
from .residue import Residue

N_PROC = get_cpu_count()
STEP_SIZE = 2000


class Chain:
    """Represents a single polypeptide chain within a PDB entry.

    Holds residues, alignment sequence, residue maps, and segment mappings
    between the PDB chain and UniProt isoforms.
    """

    def __init__(
        self,
        mmcif: mmcif_helper.mmCIF,
        pdbid: str,
        auth_asym_id: str,
        entity_id: str,
        struct_asym_id: str,
        sequence: str,
        parent: object | None,
    ):
        """Initialise a Chain from mmCIF data.

        Args:
            mmcif: Parsed mmCIF object for the parent entry.
            pdbid: PDB identifier (e.g. ``"1abc"``).
            auth_asym_id: Author chain ID as deposited in the PDB.
            entity_id: mmCIF entity identifier for this chain.
            struct_asym_id: Structural asymmetric unit chain ID.
            sequence: One-letter amino acid sequence for this chain.
            parent: Parent Entry object; ``None`` if standalone.
        """
        self.mmcif = mmcif
        self.pdbid = pdbid
        self.auth_asym_id = auth_asym_id
        self.entity_id = entity_id
        self.struct_asym_id = struct_asym_id
        self.sequence = sequence
        self.parent = parent
        self.ec = None
        self.mappings = {}
        self.best = {}
        self.residues: list[Residue] = []
        self.alignment_sequence = []
        self.is_chimera = False
        self.residue_maps = {}
        self.segments = {}
        self.canonicals = []
        self.scores = {}
        self.seg_scores = {}
        self.expression_tag_start = None
        self.connectivity_mode = False

        self.ec = mmcif.get_ec(self.entity_id)

        self.tax_fix = TaxonomyFix()
        # skips the chain
        # (e.g. MODE_NF90 and the chain is a chimera,
        # coverage < NF90_coverage, etc.)
        self.skip = False

        for k, v in mmcif.get_residues(self.auth_asym_id):
            if isinstance(v, tuple):
                tmp = Residue(
                    k, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]
                )
            else:
                tmp = []
                for v_sub in v:
                    tmp.append(
                        Residue(
                            k,
                            v_sub[0],
                            v_sub[1],
                            v_sub[2],
                            v_sub[3],
                            v_sub[4],
                            v_sub[5],
                            v_sub[6],
                            v_sub[7],
                            v_sub[8],
                        )
                    )

            self.residues.append(tmp)

        # Use the UniProt if chimera
        if self.is_chimera:
            self.tax_id = None
        # Use the mmCIF otherwise
        else:
            self.tax_id = mmcif.get_tax(self.auth_asym_id)

            # if mmCIF is None then get from taxonomy_fix
            if self.tax_id is None:
                xx = self.tax_fix.get(self.pdbid, self.entity_id)
                if xx is not None:
                    self.tax_id = xx

        self.__gen_alignment_sequence()

    # Generate the "alignment sequence": an alternative version
    # of the original sequence where several modifications are made
    # in order to simplify the alignment
    def __gen_alignment_sequence(self) -> None:
        """Build the alignment sequence from residues.

        Replaces insertions, linkers and expression tags with ``'J'``,
        uses the original residue letter for engineered mutations and
        conflicts, and expands chromophore residues into their constituent
        one-letter codes.  Result is stored in ``self.alignment_sequence``.
        No-op if the sequence has already been generated.
        """
        if self.alignment_sequence != []:
            return

        for r in self.residues:
            if isinstance(r, list):
                r = r[0]

            # Remove:
            #   initial M

            # if self.alignment_sequence == [] and r.oneL == 'M':
            #    self.alignment_sequence.append('Z')

            # "Remove" and get the begining of the expression tag to remove methionine before expression tag:
            #   Insertion
            #   Linker
            #   Expression Tag

            if (
                r.rtype == "Expression tag"
                and self.expression_tag_start is None
            ):
                self.expression_tag_start = r.n

            if r.rtype in ("Insertion", "Linker", "Expression tag"):
                self.alignment_sequence.append("J")

            # "Ignore":
            #   Engineered mutations
            #   Conflict
            #   Cloning artifact

            elif r.rtype in (
                "Engineered mutation",
                "Conflict",
                "Cloning artifact",
            ):
                self.alignment_sequence.append(
                    r.oneL_original if r.oneL_original else "X"
                )

            # unfold the chromophores
            elif r.is_chromophore:
                self.alignment_sequence.extend(list(r.oneL))

            else:
                self.alignment_sequence.append(r.oneL)

    # Compensate for index deviations due to one letter codes
    # with more than one character
    def mod_after_alignment(self, al: tuple) -> str:
        """Correct alignment positions for multi-character one-letter codes.

        Some residues (e.g. chromophores) are represented by more than one
        character in the one-letter code.  After alignment, start/stop indices
        need to be adjusted so they refer to the correct residue positions.

        Args:
            al: A two-element alignment tuple ``(unp_seq, pdb_seq)`` as
                returned by the lalign36 wrapper.

        Returns:
            The corrected PDB aligned sequence as a plain string.
        """
        seq = al[1].seq
        out = ""

        i = 0

        while i < (al[1]._al_start - 1):
            r = self.residues[i]

            if isinstance(r, list):
                r = r[0]

            if len(r.oneL) > 1:
                al[1]._al_start -= len(r.oneL) - 1
                al[1]._al_stop -= len(r.oneL) - 1
                # i -= (len(r.oneL) - 1)

            i += 1

        i = al[1]._al_start - 1
        j = 0

        while j < len(seq):
            s = seq[j]

            if s == "-":
                out += s
                j += 1
                continue

            if isinstance(self.residues[i], list):
                r = self.residues[i][0]
            else:
                r = self.residues[i]

            out += r.oneL

            i += 1
            j += len(r.oneL)

        al[1].seq = Seq(out)
        return out

    # generate the residue_maps between the chain and each isoform
    def get_each_resmap(self, row: tuple) -> tuple:
        """Compute the residue map for a single UniProt isoform.

        Iterates over the pairwise alignment for one isoform and builds a
        dictionary mapping PDB sequence positions (1-based) to UniProt
        sequence positions.  Runs connectivity refinement when enabled.

        Args:
            row: ``(isoform_accession, list_of_mappings)`` tuple taken from
                ``self.mappings.items()``.

        Returns:
            ``(isoform_accession, residue_map_dict)`` where the dict maps
            PDB position (int) → UniProt position (int or list[int] for
            chromophores).
        """
        iso, mappings = row
        residue_map = {}
        # Process in reverse so the best mappings is
        # overriding the rest
        for m in mappings[::-1]:
            al = m[2]
            if not self.is_chimera and self.connectivity_mode:
                refined_al = ConnectivityCheck(self)
                al = refined_al.alignment_refining(al, iso)
            pdb_seq = al[1]._seq
            pdb_start = al[1]._al_start
            pdb_stop = al[1]._al_stop
            unp_seq = al[0]._seq
            unp_start = al[0]._al_start
            unp_stop = al[0]._al_stop
            pdb_i = 0
            unp_i = 0
            pdb_shift = 0
            unp_shift = 0
            while (
                pdb_i + pdb_start <= pdb_stop and unp_i + unp_start <= unp_stop
            ):
                pdb_r = pdb_seq[pdb_i + pdb_shift]
                unp_r = unp_seq[unp_i + unp_shift]
                r = self.residues[pdb_i + pdb_start - 1]
                if isinstance(r, list):
                    r = r[0]
                # Remove undesired mappings just in case they were included
                if (
                    unp_r != "-"
                    and pdb_r != "-"
                    and r.rtype in ("Insertion", "Linker", "Expression tag")
                ):
                    pdb_i += 1
                    unp_i += 1
                    continue
                # Remove undesired met before expression tag to prevent single aa segment
                if (
                    r.rtype == "Initiating methionine"
                    and self.expression_tag_start is not None
                    and r.n < self.expression_tag_start
                ):
                    pdb_i += 1
                    unp_i += 1
                    continue
                if len(r.oneL) > 1 and pdb_r != "-" and unp_r != "-":
                    residue_map[pdb_i + pdb_start] = []
                    for x in range(len(r.oneL)):
                        residue_map[pdb_i + pdb_start].append(
                            unp_i + unp_start + x
                        )
                    unp_i += 1
                    pdb_i += 1
                elif pdb_r != "-" and unp_r != "-":
                    residue_map[pdb_i + pdb_start] = unp_i + unp_start
                    pdb_i += 1
                    unp_i += 1
                elif pdb_r == "-":
                    unp_i += 1
                    pdb_shift += 1
                elif unp_r == "-":
                    pdb_i += 1
                    unp_shift += 1
                if len(r.oneL) > 1:
                    pdb_shift += len(r.oneL) - 1
                    unp_i += len(r.oneL) - 1

        self.residue_maps[iso] = residue_map
        return iso, residue_map

    # generate the residue_map between the chain and each isoform
    def generate_residue_maps(self) -> None:
        """Generate residue maps for all isoforms in parallel.

        Dispatches ``get_each_resmap`` across all entries in
        ``self.mappings`` using a thread pool.  After completion, removes
        overlapping residues for chimeric chains and builds ``self.segments``
        from the resulting residue maps.
        """
        # create a process pool that uses all cpus
        with Pool(N_PROC) as pool:
            # call the function for each item in parallel, get results as tasks complete
            my_list = list(self.mappings.items())
            list(
                tqdm.tqdm(
                    pool.imap_unordered(
                        self.get_each_resmap, my_list, chunksize=STEP_SIZE
                    ),
                    disable=True,
                )
            )

        # remove the residues which map to more than one accession
        # keeping the ones that benefit continuity
        if self.is_chimera:
            self.__overlapping_residues()

        # generate the segments once we have the residue maps
        self.__segments_from_residues()

    # remove the residues which map to more than one accession
    # keeping the ones that benefit continuity
    # (only for chimeras)
    def __overlapping_residues(self) -> None:
        """Remove residues that map to multiple accessions (chimera handling).

        For each pair of isoform maps, residues present in both are resolved
        by preferring the mapping without a sequence conflict (i.e. where the
        PDB residue matches the UniProt residue).  If both match equally, the
        mapping that benefits positional continuity is kept.  Only called for
        chimeric chains.
        """
        iso_list = list(self.residue_maps.keys())

        for idx, iso1 in enumerate(iso_list):
            for iso2 in iso_list[idx + 1 :]:
                maps1 = self.residue_maps[iso1]
                maps2 = self.residue_maps[iso2]

                for key in list(maps1.keys()):
                    if key in maps2:
                        # The mapping without conflict has preference
                        pdb_r = self.sequence[key - 1]
                        try:
                            unp1_r = self.parent.accessions[iso1].seq_isoforms[
                                iso1
                            ][maps1[key] - 1]
                            unp2_r = self.parent.accessions[iso2].seq_isoforms[
                                iso2
                            ][maps2[key] - 1]

                            # If one is a conflict and the other one is not
                            if unp1_r != unp2_r and pdb_r in (unp1_r, unp2_r):
                                if pdb_r == unp2_r:
                                    del maps1[key]
                                else:
                                    del maps2[key]

                                continue
                        except TypeError:
                            # maps1[key] and/or maps2[key] are lists
                            # (the residue is a chromophore)
                            pass

                        # Otherwise try to keep continuity
                        if key - 1 in list(maps1.keys()):
                            del maps2[key]
                        else:
                            del maps1[key]

    def __group_elements(self, lst) -> list:
        """Group consecutive (PDB pos, UNP pos) pairs into segment ranges.

        Args:
            lst: Iterable of ``(pdb_position, unp_position)`` tuples, sorted
                by PDB position.

        Returns:
            List of ``((pdb_start, pdb_end), (unp_start, unp_end))`` tuples
            representing contiguous segments.
        """
        ranges = []

        for _, g in groupby(enumerate(lst), lambda i_x: i_x[0] - i_x[1][0]):
            group = list(map(itemgetter(1), g))
            ranges.append(
                ((group[0][0], group[-1][0]), (group[0][1], group[-1][1]))
            )

        # get only the min/max if there is a 1 to many mapping
        for idx, r in enumerate(ranges):
            if isinstance(r[1][0], list) or isinstance(r[1][1], list):
                start = min(r[1][0]) if isinstance(r[1][0], list) else r[1][0]
                end = max(r[1][1]) if isinstance(r[1][1], list) else r[1][1]

                ranges[idx] = (r[0], (start, end))

        return ranges

    def __overlapping(self, r1: tuple, r2: tuple) -> bool:
        """Return True if range *r1* overlaps with range *r2*.

        Args:
            r1: ``(start, end)`` tuple.
            r2: ``(start, end)`` tuple to compare against.

        Returns:
            ``True`` if the two ranges share at least one position.
        """
        return r2[0] <= r1[0] <= r2[1] or r2[0] <= r1[1] <= r2[1]

    def __segments_from_residues(self) -> None:
        """Build ``self.segments`` from residue maps and resolve chimera overlaps.

        Calls ``__group_elements`` for each isoform's residue map to produce
        contiguous segment ranges, then removes overlapping segments for
        chimeric chains so that each PDB position is assigned to at most one
        UniProt accession.
        """
        for iso, maps in list(self.residue_maps.items()):
            self.segments[iso] = self.__group_elements(
                (x, y) for x, y in sorted(maps.items(), key=itemgetter(0))
            )

        # remove overlapping segments for chimeras.
        # we create a final_set of non-overlapping segments and
        # mark all the ones overlapping with it for removal
        # by inserting them in to_remove
        if self.is_chimera:
            final_set = []
            to_remove = {}

            for iso, maps in list(self.segments.items()):
                for m, u in maps:
                    for r in final_set:
                        if self.__overlapping(m, r):
                            to_remove.setdefault(iso, []).append((m, u))
                            break

                    if m not in list(to_remove.values()):
                        final_set.append(m)

            for key, val in list(to_remove.items()):
                for seg in val:
                    self.segments[key].remove(seg)

    def get_residue_auth(self, n: int) -> tuple[str | None, str | None]:
        """Return the author (PDB) residue number and insertion code for a seq position.

        Args:
            n: 1-based sequence position in the mmCIF ``_pdbx_poly_seq_scheme``.

        Returns:
            ``(auth_seq_num, insertion_code)`` tuple.  Either value is ``None``
            when the residue is unobserved or the field is absent.  The
            insertion code defaults to a single space ``" "`` when present but
            empty.
        """
        unk = (".", "?", None, False)

        for r in self.residues:
            if isinstance(r, list):
                r = r[0]

            if r.n == n:
                return (
                    r.auth_n if r.auth_n not in unk else None,
                    r.auth_ins if r.auth_ins not in unk else " ",
                )

        return (None, None)

    def __repr__(self) -> str:
        """Return a human-readable representation of the chain."""

        return f"{self.pdbid}_{self.auth_asym_id} (E:{self.entity_id}|S:{self.struct_asym_id})"
