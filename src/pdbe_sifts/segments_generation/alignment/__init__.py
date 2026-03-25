import io
import os
import subprocess
import tempfile

from Bio import AlignIO

from pdbe_sifts.base.log import logger

GAP_OPEN = 10
GAP_EXTEND = 2

PWD = os.path.dirname(os.path.realpath(__file__))
LALIGN = "lalign36"


def do_alignment_lalign36(seq1, seq2, gap_open=GAP_OPEN, gap_extend=GAP_EXTEND):
    """Run a local pairwise alignment with lalign36 and return all alignments.

    Writes both sequences to temporary FASTA files, invokes the ``lalign36``
    binary with ``-J -m 10`` (machine-readable output), and parses the result
    with :func:`Bio.AlignIO.parse`.  Sequences shorter than 10 residues are
    padded with ``X`` to satisfy lalign36's minimum input length.

    Args:
        seq1: First sequence string (typically the UniProt sequence).
        seq2: Second sequence string or list of characters (typically the
            PDB alignment sequence).
        gap_open: Gap-opening penalty passed to lalign36 via ``-f``.
            Defaults to ``GAP_OPEN`` (10).
        gap_extend: Gap-extension penalty passed to lalign36 via ``-g``.
            Defaults to ``GAP_EXTEND`` (2).

    Returns:
        A :class:`Bio.AlignIO` iterator yielding
        :class:`Bio.Align.MultipleSeqAlignment` objects, one per local
        alignment reported by lalign36.
    """
    # cast into string if we get a list
    if isinstance(seq2, list):
        seq2 = "".join(seq2)

    # if the sequence is too short for lalign36, complete with 'X'
    if len(seq2) < 10:
        seq2 += "X" * (10 - len(seq2))

    with (
        tempfile.NamedTemporaryFile(mode="wt") as f1,
        tempfile.NamedTemporaryFile(mode="wt") as f2,
    ):
        f1.write(f">unp\n{seq1}")
        f2.write(f">pdb\n{seq2}")

        f1.flush()
        f2.flush()
        # "%s -s VT10 -J -m 10 -f %d -g %d %s %s" % (
        cmd = [
            LALIGN,
            "-J",
            "-m",
            "10",
            "-f",
            str(gap_open),
            "-g",
            str(gap_extend),
            f1.name,
            f2.name,
        ]
        logger.debug(cmd)
        output = subprocess.check_output(
            cmd, encoding="utf8", stderr=subprocess.DEVNULL
        )

        return AlignIO.parse(io.StringIO(output), "fasta-m10")


def hard_wrap(text, width=60):
    """Yield successive fixed-width chunks of *text*.

    Args:
        text: String to split into chunks.
        width: Maximum number of characters per chunk.  Defaults to 60.

    Yields:
        Non-overlapping substrings of *text*, each at most *width*
        characters long.
    """
    while text:
        yield text[:width]
        text = text[width:]


def annotate_alignment(seq1, seq2):
    """Produce a human-readable three-line alignment annotation string.

    Each aligned column is annotated as:
    - ``":"`` — identical residues
    - ``"*"`` — mismatch
    - ``"?"`` — either residue is ``"X"`` (unknown)
    - ``" "`` — gap in at least one sequence

    Output is hard-wrapped at 60 columns with ``UNP``, ``ALG``, and ``PDB``
    line labels.

    Args:
        seq1: First aligned sequence string (UniProt, gaps as ``"-"``).
        seq2: Second aligned sequence string (PDB, gaps as ``"-"``).

    Returns:
        Multi-line annotation string suitable for debug logging.
    """
    star = [" "] * len(seq1)

    out = []

    for idx, t in enumerate(zip(seq1, seq2, strict=False)):
        s1, s2 = t

        if s1 == "-" or s2 == "-":
            pass
        elif s2 == "X" or s1 == "X":
            star[idx] = "?"
        elif s1 != s2:
            star[idx] = "*"
        else:
            star[idx] = ":"

    for b, m, e in zip(
        hard_wrap(seq1), hard_wrap("".join(star)), hard_wrap(seq2), strict=False
    ):
        out.append(f"UNP: {b}\nALG: {m}\nPDB: {e}\n")

    return "\n".join(out)


def remove_range_alignment(pdb_ranges, unp_ranges, rm, feature, wrong_align):
    """Remove a sub-range from a set of PDB/UniProt alignment range pairs.

    Trims or splits existing range pairs so that the region *rm* is excluded
    from the PDB side.  For non-insertion features the corresponding UniProt
    positions are also removed; for insertions the UniProt numbering is
    preserved from the start of the gap.

    Args:
        pdb_ranges: List of ``(start, end)`` integer tuples covering PDB
            sequence positions.
        unp_ranges: List of ``(start, end)`` integer tuples covering the
            corresponding UniProt positions (same length as *pdb_ranges*).
        rm: ``(start, end)`` tuple identifying the PDB sub-range to remove.
        feature: Annotation type string (e.g. ``"Insertion"``); controls
            whether UniProt positions are renumbered after the gap.
        wrong_align: Reserved flag (currently unused in the range arithmetic).

    Returns:
        ``(new_pdb_ranges, new_unp_ranges)`` tuple with the excised region
        excluded.
    """
    new_pdb = []
    new_unp = []

    remove = feature != "Insertion"  # or wrong_align
    for pdb, unp in zip(pdb_ranges, unp_ranges, strict=False):
        if rm[0] >= pdb[0]:
            if rm[0] > pdb[1]:
                new_pdb.append(pdb)
                new_unp.append(unp)
                continue
            if rm[0] > pdb[0]:
                new_pdb.append((pdb[0], rm[0] - 1))
                new_unp.append((unp[0], unp[0] + rm[0] - pdb[0] - 1))
            if rm[1] < pdb[1]:
                new_pdb.append((rm[1] + 1, pdb[1]))
                if remove:
                    new_unp.append((unp[0] + rm[1] - pdb[0] + 1, unp[1]))
                else:
                    new_unp.append((unp[0] + rm[0] - pdb[0], unp[1]))

        else:
            if rm[1] < pdb[0]:
                new_pdb.append(pdb)
                new_unp.append(unp)
                continue
            if pdb[1] > rm[1] > pdb[0]:
                new_pdb.append((rm[1] + 1, pdb[1]))
                if remove:
                    new_unp.append((unp[0] + rm[1] - pdb[0] + 1, unp[1]))
                else:
                    new_unp.append((unp[0] + rm[0] - pdb[0], unp[1]))

    return (new_pdb, new_unp)


def get_conflicts(seq1, seq2):
    """Count mismatched positions between two aligned sequences.

    Args:
        seq1: First aligned sequence string.
        seq2: Second aligned sequence string.

    Returns:
        Number of positions where the characters differ (gap characters are
        compared as literal characters).
    """
    return sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2, strict=False))


def get_identity(seq1, seq2):
    """Compute the fractional sequence identity of an alignment.

    Identity is calculated as the number of identical positions divided by
    the length of *seq2* (the PDB sequence), rounded to two decimal places.

    Args:
        seq1: First aligned sequence string (UniProt).
        seq2: Second aligned sequence string (PDB); its length is used as
            the denominator.

    Returns:
        Fractional identity in the range ``[0.0, 1.0]``, rounded to 2 d.p.
    """
    return round(
        sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2, strict=False))
        / float(len(seq2)),
        2,
    )


def get_score(seq1, seq2):
    """Return the alignment score for a sequence pair.

    Currently delegates directly to :func:`get_identity`, so the score is
    the fractional sequence identity normalised by the PDB sequence length.

    Args:
        seq1: First aligned sequence string (UniProt).
        seq2: Second aligned sequence string (PDB).

    Returns:
        Fractional identity score in the range ``[0.0, 1.0]``.
    """
    return get_identity(seq1, seq2)


def get_align_chunk(seq, start, stop):
    """Extract a sub-sequence from an aligned string using non-gap indexing.

    Iterates over *seq*, ignoring gap characters (``"-"``), and returns the
    slice corresponding to 1-based positions *start* through *stop* of the
    ungapped sequence together with any gap characters interspersed between
    those positions.

    Args:
        seq: Aligned sequence string, possibly containing ``"-"`` gap
            characters.
        start: 1-based start position in the ungapped sequence (inclusive).
        stop: 1-based stop position in the ungapped sequence (inclusive).

    Returns:
        Substring of *seq* covering positions *start*–*stop*, including any
        internal gap characters.  Returns an empty string when
        ``start > stop``.
    """
    if start > stop:
        return ""

    out = ""
    i = 0

    for c in seq:
        if c != "-":
            i += 1
        if i >= start:
            out += c
        if i == stop:
            return out

    return out
