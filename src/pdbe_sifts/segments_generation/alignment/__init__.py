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
    while text:
        yield text[:width]
        text = text[width:]


def annotate_alignment(seq1, seq2):
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
    return sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2, strict=False))


def get_identity(seq1, seq2):
    return round(
        sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2, strict=False))
        / float(len(seq2)),
        2,
    )


def get_score(seq1, seq2):
    return get_identity(seq1, seq2)


# Get a chunk of the sequence from start to stop (1-indexed)
# ignoring annotation gaps (i.e. '-').
def get_align_chunk(seq, start, stop):
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
