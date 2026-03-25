import contextlib
from enum import Enum
from multiprocessing.dummy import Pool
from typing import NamedTuple

import tqdm
from funcy.debug import log_durations

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import get_cpu_count
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.segments_generation import alignment

# from segments_gen.uniref90_pkl import NF90Coverage, NF90TaxID
from pdbe_sifts.unp.unp import UNP, get_unp_object

NF_COVERAGE = 0.7
N_PROC = get_cpu_count()
STEP_SIZE = 2000


class SMapping(NamedTuple):
    """A single UniProt segment mapping for one PDB chain.

    Attributes:
        accession: UniProt accession string (e.g. ``"P12345"``).
        range_start: First UniProt residue position covered by this mapping
            (1-based, inclusive).
        range_stop: Last UniProt residue position covered by this mapping
            (1-based, inclusive).
    """

    accession: str
    range_start: int
    range_stop: int


class CustomSequenceAccession:
    """Minimal UNP-like object for user-provided FASTA sequences.

    Allows arbitrary protein sequences to be used in the SIFTS alignment
    pipeline without requiring a UniProt accession.  The interface mirrors
    the subset of UNP attributes that are actually accessed by helper.py
    and generate_xref_csv.py.
    """

    def __init__(self, accession: str, sequence: str, name: str = ""):
        """Initialise a custom sequence accession object.

        Args:
            accession: Identifier string used in FASTA headers and as the
                isoform key (e.g. a chain ID or arbitrary label).
            sequence: One-letter amino-acid sequence.
            name: Human-readable display name.  Defaults to *accession* when
                empty.
        """
        self.accession = accession
        self.ad_dbref_auto_acc = accession  # used in get_curated_db_mappings
        self.sequence = sequence
        self.seq_isoforms = {accession: sequence}
        self.taxonomy = []  # no taxonomy info for custom sequences
        self.longName = name or accession  # fallback: use accession as name
        self.date_seq_update = [
            "",
            "0",
        ]  # [modified_date, version] — no versioning for custom seqs

    def getAllIsoforms(self) -> dict[str, str]:  # noqa: N802 — matches UNP method name
        """Return a dict of all isoforms for this custom sequence.

        For user-provided sequences there is only a single "isoform" — the
        sequence itself — keyed by the accession.  This method exists to
        satisfy the interface expected by the alignment pipeline, which calls
        ``getAllIsoforms()`` on both :class:`UNP` and custom sequence objects.

        Returns:
            Dict mapping ``{accession: sequence}``.
        """
        return {self.accession: self.sequence}


class SiftsMode(Enum):
    """Operating mode for the SIFTS alignment pipeline.

    Attributes:
        ISOFORM: Align against all UniProt isoforms of each accession.
        NF90: Align against UniRef90 representative sequences instead of
            isoforms, used for large-scale / reduced-redundancy runs.
    """

    ISOFORM = 1
    NF90 = 2


def get_accession(entry: Entry, acc: str) -> UNP:
    """Fetch or cache a UNP object for a given accession.

    First checks the entry-level accession cache.  If not present, resolves
    the accession via :func:`~pdbe_sifts.unp.unp.get_unp_object` and stores
    the result in the cache for subsequent calls.

    Args:
        entry: :class:`~pdbe_sifts.mmcif.entry.Entry` whose ``accessions``
            dict acts as the cache.
        acc: UniProt accession string to look up.

    Returns:
        A :class:`~pdbe_sifts.unp.unp.UNP` object, or ``None`` if the
        accession cannot be resolved.
    """
    if acc in entry.accessions:
        return entry.accessions[acc]

    unp = get_unp_object(acc)
    if unp:
        entry.accessions[unp.accession] = unp

    return unp


def check_range(real_ranges: list, pdb: list, unp: list) -> bool:
    """Verify that every (pdb_pos, unp_pos) pair is present in real_ranges.

    Args:
        real_ranges: List of ``(pdb_pos, unp_pos)`` pairs considered valid.
        pdb: Sequence of PDB positions to test.
        unp: Sequence of UniProt positions to test (same length as *pdb*).

    Returns:
        ``True`` if every zipped ``(pdb, unp)`` pair appears in *real_ranges*,
        ``False`` otherwise.
    """
    return all((p, u) in real_ranges for p, u in zip(pdb, unp, strict=False))


def fmt_ranges(ranges: list) -> list[list[int]]:
    """Normalise a list of ranges to ``[[start, end], ...]`` integer lists.

    Accepts ranges either as ``"start-end"`` strings (e.g. from AD_DBREF /
    Redis) or as existing integer sequences and converts all of them to
    two-element ``[int, int]`` lists.

    Args:
        ranges: List of ranges, each either a ``"start-end"`` string or an
            iterable of two integers.

    Returns:
        List of ``[start, end]`` integer lists.
    """
    out = []

    for r in ranges:
        if isinstance(r, str):
            out.append([int(x) for x in r.split("-")])
        else:
            out.append(r)

    return out


def overlapping(ranges: list) -> bool:
    """Return True if any two ranges in the list overlap.

    Normalises input via :func:`fmt_ranges` before comparison.  Two ranges
    overlap when either endpoint of the first falls within the second.

    Args:
        ranges: List of ranges in any format accepted by :func:`fmt_ranges`.

    Returns:
        ``True`` if at least one pair of ranges overlaps, ``False`` otherwise.
    """
    ranges = fmt_ranges(ranges)

    for idx, r1 in enumerate(ranges):
        for r2 in ranges[idx + 1 :]:
            if r2[0] <= r1[0] <= r2[1] or r2[0] <= r1[1] <= r2[1]:
                return True

    return False


def out_of_order(ranges: list) -> bool:
    """Return True if any range starts at or before a preceding range.

    Detects when segments are not in strictly ascending order, which
    indicates a repeated or rearranged mapping that must be handled
    specially during alignment.  Input is normalised via :func:`fmt_ranges`.

    Args:
        ranges: List of ranges in any format accepted by :func:`fmt_ranges`.

    Returns:
        ``True`` if any later range has a start position less than or equal
        to the start of an earlier range, ``False`` otherwise.
    """
    ranges = fmt_ranges(ranges)

    for idx, r1 in enumerate(ranges):
        for r2 in ranges[idx + 1 :]:
            if r2[0] <= r1[0]:
                return True

    return False


def unroll_map(m: list) -> tuple:
    """Unpack a serialised mapping record into its constituent parts.

    Mapping records may arrive as a single-element list (chain only) or as a
    three-element list ``[chain, accs_str, ranges_str]`` where the second and
    third elements are ``eval``-able string representations of Python objects.

    Args:
        m: List with either one element ``[chain]`` or three elements
            ``[chain, accs_str, ranges_str]``.

    Returns:
        ``(chain, accs, ranges)`` tuple where *accs* and *ranges* are the
        deserialised Python objects, or ``(chain, None, None)`` for a
        single-element input.
    """
    chain = m[0]

    if len(m) > 1:
        accs = eval(m[1])
        ranges = eval(m[2])
    else:
        accs = None
        ranges = None

    return chain, accs, ranges


# Currently:
#   GAP_OPEN = 10
#   GAP_EXTEND = 2
#
# Therefore, gap won't open if (AT LEAST!):
# G > E / 2 - 5
#
# G = gap length
# E = extra aligned sequence
#
# This method identifies manually annotated large gaps so we can treat them
# as overlapping mappings (forcing extra alignments, one per range)
def large_gap(ranges: list) -> bool:
    """Detect whether the gap between two consecutive ranges is unusually large.

    Uses the lalign36 gap scoring heuristic: with ``GAP_OPEN=10`` and
    ``GAP_EXTEND=2``, a gap of length *G* will only be opened if the extra
    aligned sequence *E* satisfies ``G > E / 2 - 5``.  Here a factor of 3 is
    applied to be conservative.  Ranges may be provided as ``"start-end"``
    strings (from AD_DBREF) or as integer lists.

    Only the first pair of ranges (index 0 and 1) is evaluated.

    Args:
        ranges: List of at least two ranges, each either a ``"start-end"``
            string or a two-element integer list/tuple.

    Returns:
        ``True`` if the gap between the first two ranges exceeds the
        heuristic threshold, ``False`` if fewer than two ranges are given or
        the gap is within the expected bounds.
    """
    if len(ranges) < 2:
        return False

    # If ranges come from mmCIF: array of ints
    # If ranges come from AD_DBREF (redis): array of strings
    if isinstance(ranges[0], str):
        iranges = [list(map(int, x.split("-"))) for x in ranges]
    else:
        iranges = ranges

    gap = iranges[1][0] - iranges[0][1] - 1
    extra = iranges[1][1] - iranges[1][0] + 1

    return gap > 3 * (extra / 2 - 5)


class EntryMapping:
    """Orchestrate the mapping between a single PDB chain and UniProt accessions.

    For each candidate accession, runs lalign36 alignments against every
    UniProt isoform, stores the best mapping, builds residue maps and
    segment ranges via the parent :class:`Chain` object, and optionally
    applies connectivity correction.
    """

    def __init__(
        self,
        entry: Entry,
        chain: str,
        chain_mapping: list[SMapping],
        nf90_mode: bool = False,
        NFT: dict | None = None,
        NFC: dict | None = None,
        connectivity_mode: bool = True,
    ) -> None:
        """Initialise the mapping for one chain.

        Args:
            entry: Parsed :class:`Entry` object for the PDB entry.
            chain: Author chain ID to process.
            chain_mapping: Pre-computed list of :class:`SMapping` objects
                (accession + range) for this chain, e.g. from DuckDB.
            nf90_mode: When ``True`` uses UniRef90 representative sequences
                instead of all isoforms.
            NFT: NF90 taxonomy lookup object (required when *nf90_mode* is
                ``True``).
            NFC: NF90 coverage lookup object (required when *nf90_mode* is
                ``True``).
            connectivity_mode: When ``True`` activates the connectivity
                correction step on the chain.
        """
        self.entry = entry
        self.chain = chain
        self.chain_mapping = chain_mapping
        self.nf90_mode = nf90_mode
        self.chain_obj = self.entry.chains[self.chain]
        self.repeated_acc = False
        self.accs: list[str] = []
        self.ranges: list[tuple[int, int]] = []
        self.NFT = NFT
        self.NFC = NFC
        self.connectivity_mode = connectivity_mode

    def _get_accessions(self) -> None:
        """Populate ``self.accs`` and ``self.ranges`` from mappings or mmCIF.

        First validates the provided *chain_mapping* accessions against
        UniProt.  If none are valid, falls back to accessions embedded in
        the mmCIF ``_struct_ref`` category.  Ranges are similarly taken
        from the mmCIF when not already available.
        """
        # Validate the provided mappings
        for mapping in self.chain_mapping:
            unp = get_accession(self.entry, mapping.accession)
            if unp:
                logger.info(
                    f"Valid mapping: {unp.accession} {mapping.range_start}-{mapping.range_stop}"
                )
                self.accs.append(unp.accession)
                self.ranges.append((mapping.range_start, mapping.range_stop))

        # If no valid mappings were provided, try to get them from the mmCIF
        if not self.accs:
            cif_accs = self.entry.mmcif.get_unp(self.chain)
            for acc in cif_accs:
                unp = get_accession(self.entry, acc)
                if unp:
                    self.accs.append(unp.accession)

            if self.accs:
                logger.info(f"Accessions come from mmCIF: {self.accs}")

        # If we have mappings but no ranges, try to get them from the mmCIF
        if self.accs and not self.ranges:
            self.ranges = self.entry.mmcif.get_ranges(self.chain, self.accs[0])
            logger.info(f"Ranges come from mmCIF: {self.ranges}")

    def set_chain_accessions(self) -> bool:
        """Resolve accessions and verify the chain is ready for alignment.

        Calls ``_get_accessions()``, then checks that the chain object
        exists and, in NF90 mode, that UniProt coverage meets the threshold.

        Returns:
            ``True`` if the chain can be processed, ``False`` if it should
            be skipped (no valid accession, unknown chain, or insufficient
            NF90 coverage).
        """
        self._get_accessions()
        if not self.accs:
            logger.warning("The mmCIF doesnt have a valid UniProt accession")
            return False

        logger.debug(f"Accessions [{self.chain}]: {self.accs}")
        try:
            self.chain_obj = self.entry.chains[self.chain]
        except KeyError:
            logger.warning(
                f"The chain {self.chain} was not found in the mmCIF  or it is not a polypeptide"
            )
            return False

        return not (self.nf90_mode and not self.check_unp_coverage())

    @log_durations(logger.debug)
    def process(self) -> None:
        """Run the full alignment and segment generation pipeline for this chain.

        For each unique accession:

        1. Fetches the UniProt object and all its isoforms (or NF90 reps).
        2. Aligns each isoform against the chain alignment sequence in
           parallel using lalign36.
        3. Stores alignment results and updates the best mapping.

        After processing all accessions, triggers ``generate_residue_maps``
        on the chain object to produce final segment data.
        """
        self.seq_pdb = self.chain_obj.alignment_sequence
        self.check_repeated_accession()

        logger.info(f"Processing {self.chain}: {self.accs}")

        for acc in set(self.accs):
            self.acc = acc
            isoforms = {}
            unp = get_accession(self.entry, acc)
            if not unp:
                continue
            self.unp = unp
            # isoform, score, length
            self.chain_obj.best[unp.accession] = (None, 0.0, 0)

            # Add the canonical
            self.chain_obj.canonicals.append(unp.accession)

            # Don't store isoforms if the chain is a chimera
            if self.chain_obj.is_chimera:
                isoforms[unp.accession] = unp.sequence
            elif self.nf90_mode:
                isoforms = self.get_nf90_isoforms(unp, isoforms)
            else:
                for iso, seq in list(unp.getAllIsoforms().items()):
                    isoforms[iso] = seq
                    self.entry.accessions[iso] = unp

            # create a process pool that uses all cpus
            with Pool(N_PROC) as pool:
                # call the function for each item in parallel, get results as tasks complete
                my_list = list(isoforms.items())
                list(
                    tqdm.tqdm(
                        pool.imap_unordered(
                            self.process_each_isoform,
                            my_list,
                            chunksize=STEP_SIZE,
                        ),
                        disable=True,
                    )
                )
        if self.connectivity_mode:
            self.chain_obj.connectivity_mode = True
        self.chain_obj.generate_residue_maps()
        logger.info("Segments generated")

    def process_each_isoform(self, row: tuple) -> str:
        """Align a single isoform against the chain sequence and record results.

        Designed to run inside a thread pool via ``pool.imap_unordered``.

        Args:
            row: ``(isoform_accession, isoform_sequence)`` tuple.

        Returns:
            The isoform accession string (for tqdm progress tracking).
        """
        iso, seq = row
        # Don't repeat the canonical
        if (
            not self.nf90_mode
            and seq == self.unp.sequence
            and iso != self.unp.accession
        ):
            logger.info(f"{iso} is the canonical. Skipping...")
            return iso

        all_alns = alignment.do_alignment_lalign36(seq, self.seq_pdb)
        alns = []
        # We want as many alignments as times the accession is repeated.
        # Alignments sorted from best to worst
        if self.repeated_acc:
            for _ in range(self.accs.count(self.acc)):
                try:
                    alns.append(next(all_alns))
                except StopIteration:
                    break
        else:
            # the alignment might be impossible
            # (e.g. an isoform which doesn't have the relevant sequence fragment)
            with contextlib.suppress(StopIteration):
                alns.append(next(all_alns))
        self.process_alignments(self.unp, iso, alns)
        return iso

    def process_alignments(self, unp: UNP, iso: str, alns: list) -> None:
        """Store alignment results and update scores for one isoform.

        Args:
            unp: UniProt object for the canonical accession.
            iso: Isoform accession being processed.
            alns: List of alignment tuples returned by lalign36 (one per
                repeat of the accession in chimera/repeated-accession mode).
        """
        # Loop through as many alignments as accessions in the list
        # if it is only one accession which repeats
        for al in alns:
            pdb_start, pdb_end = (al[1]._al_start, al[1]._al_stop)
            unp_start, unp_end = (al[0]._al_start, al[0]._al_stop)

            self.chain_obj.mod_after_alignment(al)

            logger.debug(alignment.annotate_alignment(al[0].seq, al[1].seq))

            logger.debug(f"PDB: [{pdb_start}-{pdb_end}]")
            logger.debug(f"UNP: [{unp_start}-{unp_end}]")

            pdb_ranges = [(pdb_start, pdb_end)]
            unp_ranges = [(unp_start, unp_end)]

            self.chain_obj.mappings.setdefault(iso, []).append(
                (pdb_ranges, unp_ranges, al)
            )
            identity = alignment.get_identity(al[0]._seq, al[1]._seq)
            score = alignment.get_score(al[0]._seq, al[1]._seq)

            self.chain_obj.scores[iso] = (identity, score)

            self.chain_obj.seg_scores.setdefault(iso, {})[str(pdb_ranges)] = (
                identity
            )

            self._update_best_mapping(unp, iso, al, score)

    def get_nf90_isoforms(self, unp: UNP, isoforms: dict) -> dict:
        """Fetch NF90 representative sequences for a UniProt accession.

        Args:
            unp: UniProt object for the accession being processed.
            isoforms: Existing isoform dict to extend in-place.

        Returns:
            Updated *isoforms* dict with NF90 representative sequences added.
        """
        logger.debug(f"fetching NF90 isoforms for {unp}")
        for acc in self.NFT.get_nf90(unp.accession):
            unp_nf90 = get_accession(self.entry, acc)

            if not unp_nf90:
                continue

            acc = unp_nf90.accession if "-" not in acc else acc

            if acc not in unp_nf90.seq_isoforms:
                unp_nf90.seq_isoforms[acc] = unp_nf90.sequence

            if acc not in self.entry.accessions:
                self.entry.accessions[acc] = unp_nf90

            isoforms[acc] = unp_nf90.sequence
        return isoforms

    def get_valid_accession(self, acc: str) -> UNP | None:
        """Fetch and validate a UniProt accession, falling back to mmCIF data.

        If *acc* cannot be resolved directly, tries the accession embedded
        in the mmCIF ``_struct_ref`` category.  Resets chimera and
        repeated-accession flags if resolution fails entirely.

        Args:
            acc: UniProt accession to validate.

        Returns:
            A :class:`UNP` object, or ``None`` if no valid accession
            can be found.
        """
        unp = get_accession(self.entry, acc)
        if not unp:
            accs = self.entry.mmcif.get_unp(self.chain)
            if not accs:
                logger.warning(
                    f"The mmCIF doesnt have a UniProt accession: {self.entry.pdbid}"
                )
                self.chain_obj.is_chimera = False
                self.repeated_acc = False
                return

            if len(accs) > 1:
                logger.error(
                    "It is a chimera and at least one of the accessions is not valid anymore."
                )
                self.chain_obj.is_chimera = False
                self.repeated_acc = False
                return

            logger.warning(
                f"The accession {acc} is not valid. Got {accs[0]} from mmCIF"
            )
            if acc == accs[0]:
                logger.warning(
                    f"CIF has same invalid obsolete accession {acc}. Doing nothing"
                )
                return

            unp = get_accession(self.entry, accs[0])

        return unp

    def _update_best_mapping(
        self, unp: UNP, iso: str, al: tuple, score: float
    ) -> None:
        """Update the best mapping for an accession if this alignment is better.

        The best mapping is replaced when the new *score* is strictly higher,
        or equal with a longer aligned sequence, or equal/same-length but the
        isoform is canonical or matches the AD_DBREF accession.

        Args:
            unp: UniProt object for the accession.
            iso: Isoform accession being evaluated.
            al: Alignment tuple ``(unp_seq, pdb_seq)``.
            score: Alignment score for this isoform.
        """
        try:
            ad_dbref_acc = self.chain_mapping[0].accession
        except Exception:
            logger.warning(
                f"Could not get the AD_DBREF accession for {self.entry.pdbid} {self.chain}"
            )
            ad_dbref_acc = None
        best_score = self.chain_obj.best[unp.accession][1]
        best_seq_len = self.chain_obj.best[unp.accession][2]
        if (
            score > best_score
            or (score == best_score and len(al[1]._seq) > best_seq_len)
            or (score == best_score and len(al[1]._seq) == best_seq_len)
            and (iso in self.chain_obj.canonicals or iso == ad_dbref_acc)
        ):
            self.chain_obj.best[unp.accession] = (
                iso,
                score,
                len(al[1]._seq),
            )

    def check_unp_coverage(self) -> bool:
        """Check whether the chain meets the NF90 coverage threshold.

        Returns:
            ``False`` if the UniProt object is unavailable, the chain is
            chimeric or has a repeated accession, or coverage is below
            ``NF_COVERAGE`` (0.7).  ``True`` otherwise.
        """
        logger.debug("Checking coverage")
        unp = get_accession(self.entry, self.accs[0])

        if not unp:
            return False

        if self.repeated_acc or self.chain_obj.is_chimera:
            logger.warning("It is a chimera so we skip it for UniRef90")
            return False
        coverage = self.NFC.get_coverage(
            self.entry.pdbid, self.chain, unp.accession
        )
        if coverage < NF_COVERAGE:
            logger.warning(
                f"Coverage {coverage} < {NF_COVERAGE} so we skip it for UniRef90"
            )
            return False
        return True

    def check_repeated_accession(self) -> None:
        """Detect whether the same accession appears multiple times (chimera / repeat).

        Sets ``self.repeated_acc = True`` when the same accession maps to
        out-of-order, overlapping, or large-gap ranges.  Sets
        ``chain_obj.is_chimera = True`` when multiple distinct accessions
        are present for the chain.
        """
        if len(self.accs) > 1:
            # One accession that repeats
            if len(set(self.accs)) == 1:
                if out_of_order(self.ranges) or overlapping(self.ranges):
                    self.repeated_acc = True
                    logger.warning("Same accession repeats!")
                # Large gap go hereƒ
                elif large_gap(self.ranges):
                    # Calculate gap and extra aligned sequence and,
                    # if large enough, treat as a repeated_acc
                    self.repeated_acc = True
                    logger.warning("Large gap detected")
            else:
                logger.info(f"Chain {self.chain} is Chimera!")
                self.chain_obj.is_chimera = True
