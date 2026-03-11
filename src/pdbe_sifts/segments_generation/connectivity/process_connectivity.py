#!/usr/bin/env python3

import math
import regex as re

from Bio.Seq import Seq

from pdbe_sifts.segments_generation.connectivity.ccd_parser import CcdFile
from pdbe_sifts.base.log import logger
from pdbe_sifts.config import load_config

conf = load_config()

PEPTIDE_BOND_LEN = 1.42
GAP_ALIGNMENT_PATTERN = r"(?<=-)([A-Z]{1,5})(?=-)"
CONS_ALIGNMENT_PATTERN = r"[A-Za-z]+(?:-[A-Za-z]+)*"
STARTING_GAP = r"^((?:-?[A-Z]){1,5})(-{2,})(.*?)$"
ENDING_GAP = "^(.*?)(-{2,})((?:-?[A-Z]){1,5})$"
RESIDUE_PATTERN = r"([A-Z])"


def fmt_ranges(ranges):
    out = []
    for r in ranges:
        if isinstance(r, str):
            out.append([int(x) for x in r.split("-")])
        else:
            out.append(r)
    return out


def overlapping(ranges):
    ranges = fmt_ranges(ranges)
    for idx, r1 in enumerate(ranges):
        for r2 in ranges[idx + 1 :]:
            if (
                r2[0] <= r1[0] <= r2[1]
                or r2[0] <= r1[1] <= r2[1]
                or r1[0] <= r2[0] <= r1[1]
                or r1[0] <= r2[1] <= r1[1]
            ):
                return True
    return False


def merge_segment(seg1, seg2):
    start_pdb = seg1[0][0]
    end_pdb = seg2[0][1]
    start_unp = seg1[1][0]
    end_unp = seg2[1][1]
    return ((start_pdb, end_pdb), (start_unp, end_unp))


def compute_atom_distance(atom1, atom2):
    """Compute the Euclidean distance between two 3D coordinates and return a truncated value.
    The input should be a list of coordinates in the format [x, y, z]."""
    x1 = atom1[0]
    y1 = atom1[1]
    z1 = atom1[2]
    x2 = atom2[0]
    y2 = atom2[1]
    z2 = atom2[2]
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    truncated_number = int(dist * 100) / 100
    return truncated_number


def is_valid_residue_status(r1, r2):
    """Check the status of two residues and ensure the connectivity check is valid.
    Takes two residue objects as input."""
    # we don't need to process insertion, linker etc...
    to_prevent = ["Insertion", "Linker", "Chromophore"]
    if (
        r1.rtype in to_prevent
        or r2.rtype in to_prevent
        or r1.observed is False
        or r2.observed is False
        or r1.auth_ins is not False
        or r2.auth_ins is not False
    ):
        return False
    return True


class ConnectivityCheck:
    """Class responsible for performing connectivity checks at various stages of the segment generation process."""

    def __init__(self, chain_obj, repeated_acc=None):
        self.chain_obj = chain_obj
        self.repeated_acc = repeated_acc
        self.r1_ccd_ters = None
        self.r2_ccd_ters = None
        self.r1_ters_coordinates = None
        self.r2_ters_coordinates = None

    def get_residues_ccd_ters(self, r1, r2):
        """
        Function to retrieve the C-terminal and N-terminal atoms of two residues.

        Args:
            r1 (int): The index of the first residue (1 to N).
            r2 (int): The index of the second residue (1 to N).
        """
        r1 = self.chain_obj.residues[r1 - 1]
        r2 = self.chain_obj.residues[r2 - 1]
        if is_valid_residue_status(r1, r2):
            self.r1_ccd_ters = CcdFile(r1.threeL).process()
            self.r2_ccd_ters = CcdFile(r2.threeL).process()
        else:
            self.r1_ccd_ters = 0
            self.r2_ccd_ters = 0

    def get_residue_ters_coordinates(self, r1, r2):
        """
        Function to retrieve the C-terminal and N-terminal atoms 3D coordinates of two residues.

        Args:
            r1 (int): The index of the first residue (1 to N).
            r2 (int): The index of the second residue (1 to N).
        """
        try:
            r1after = self.chain_obj.residues[r1]
            r2before = self.chain_obj.residues[r2 - 2]
        except IndexError:
            r1after = self.chain_obj.residues[r1 - 1]
            r2before = self.chain_obj.residues[r2 - 1]
        r1 = self.chain_obj.residues[r1 - 1]
        r2 = self.chain_obj.residues[r2 - 1]
        if (
            self.r1_ccd_ters == 0
            or self.r2_ccd_ters == 0
            or not is_valid_residue_status(r1, r2)
            or not is_valid_residue_status(r1after, r2before)
        ):
            self.r1_ters_coordinates = None
            self.r2_ters_coordinates = None
        else:
            self.r1_ters_coordinates = (
                self.chain_obj.mmcif.get_ters_coordinates_of_residue(
                    self.chain_obj.struct_asym_id, r1.n, self.r1_ccd_ters
                )
            )
            self.r2_ters_coordinates = (
                self.chain_obj.mmcif.get_ters_coordinates_of_residue(
                    self.chain_obj.struct_asym_id, r2.n, self.r2_ccd_ters
                )
            )

    def check_res_conn(self, r1, r2):
        """
        Function to check the connectivity of two residues.

        Args:
            r1 (int): The index of the first residue (1 to N).
            r2 (int): The index of the second residue (1 to N).
        """
        # When a chromophore is present, the length of self.chain_obj.residues may differ from that of the alignment sequence.
        # This discrepancy can cause index out-of-range errors.
        # For example, in entry 3EVP, the alignment boundary is at 245, while chain.residues has only 243 elements.
        r1_shift = sum(
            [
                len(r.oneL) - 1
                for r in self.chain_obj.residues[:r1]
                if r.rtype == "Chromophore"
            ]
        )
        r2_shift = sum(
            [
                len(r.oneL) - 1
                for r in self.chain_obj.residues[:r2]
                if r.rtype == "Chromophore"
            ]
        )
        r1 = r1 - r1_shift
        r2 = r2 - r2_shift
        try:
            self.get_residues_ccd_ters(r1, r2)
            self.get_residue_ters_coordinates(r1, r2)
            return self.compute_terminals_distance()
        except IndexError:
            pass

    def boundaries_check(self, subseq, mode, r1, aligned_res_bound):
        if mode == "start":
            shift = subseq.count("-")
            limit = r1 + aligned_res_bound - shift - 1
            residues_to_check = [r.n for r in self.chain_obj.residues[: limit + 1]]
            for i in range(len(residues_to_check) - 1):
                connected = self.check_res_conn(
                    residues_to_check[i], residues_to_check[i + 1]
                )
                if not connected:
                    return subseq
            subseq = "".join(subseq)
            dash_to_add = subseq.count("-")
            residues = "-" * dash_to_add + "".join(subseq.split("-"))
        else:
            from_res = aligned_res_bound - sum([1 for r in subseq[:] if r != "-"]) - 1
            residues_to_check = [
                r.n for r in self.chain_obj.residues[from_res : aligned_res_bound + 1]
            ]
            for i in range(len(residues_to_check) - 1):
                connected = self.check_res_conn(
                    residues_to_check[i], residues_to_check[i + 1]
                )
                if not connected:
                    return subseq
            subseq = "".join(subseq)
            residues = "".join(subseq.split("-"))
        return list(residues)

    def subseq_find_midgap(
        self, pdb_seq, new_pdb_seq, res_first, res_second, start_res
    ):
        """
        Function to refine an extended gap between two continuous regions, where more than 5 PDB residues are aligned (to a UniProt residue).

        This function aims to move residues to the right continuous region based on connectivity, improving the alignment and continuity between the two regions.

        Args:
            pdb_seq (string): the aligned sequence. It contains residues and dashes.
            new_pdb_seq (list): same that pdb_seq but this is a list.
            res_first (int): residue index where the extended gap starts.
            res_second (int): residue index where the extended gap ends.
            start_res (int): residue chain number where the aligned sequence (pdb_seq) starts. Indeed, the aligned sequence doesn't always start at 0 or 1 where the chain does.

        Returns:
            new_pdb_seq: the updated sequence as a list of residues and dashes.
        """
        subseq = pdb_seq[res_first + 1 : res_second]
        # Identify residues within the subsequence and check their connectivity
        one_res_matches = [
            (match.group(), match.start(), match.end() - 1)
            for match in re.finditer(RESIDUE_PATTERN, subseq)
        ]
        if not one_res_matches:
            return new_pdb_seq
        # Starting index for alignment refinement, can be adjusted based on connectivity
        rconn_beg = 0
        # Determine the position of the first and last matched residues within the subsequence
        pos_start_subseq = one_res_matches[0][1]
        pos_end_subseq = one_res_matches[-1][1]
        # Append the next residue in the continuous region to verify connectivity with the last residue in the gap
        one_res_matches.append((pdb_seq[res_second + 1], len(subseq), len(subseq)))
        # Convert subsequence positions to full sequence positions (pdb_seq indexes i.e the aligned sequence)
        pos_start_pdbseq = res_first + pos_start_subseq + 1
        pos_end_pdbseq = res_first + pos_end_subseq + 1
        # Compute shift due to gaps ('-') to map positions to residue chain indexes
        current_shift = pdb_seq[:pos_start_pdbseq].count("-")
        last_shift = pdb_seq[:pos_end_pdbseq].count("-")
        # Convert pdb_seq positions to chain sequence positions (which includes non-aligned regions)
        pos_start_chain = pos_start_pdbseq + start_res - current_shift - 1
        pos_end_chain = pos_end_pdbseq + start_res - last_shift
        # Generate the list of residues to check for connectivity; +1 to check the residue located the first position of the next continuous region
        residues_to_check = list(range(pos_start_chain, pos_end_chain + 1))
        # Iterate backwards through residues to verify connectivity from left to right
        for i in range(len(residues_to_check) - 1, 0, -1):
            r2 = self.chain_obj.residues[residues_to_check[i]].n
            r1 = self.chain_obj.residues[residues_to_check[i - 1]].n
            connected = self.check_res_conn(r1, r2)
            if not connected:
                # If connectivity is broken, adjust alignment refinement boundaries
                pos_start_subseq = one_res_matches[i][1]
                # Adjust starting index to exclude non-connected residues
                rconn_beg = pos_start_subseq - subseq[: one_res_matches[i][1]].count(
                    "-"
                )
                # Update position in pdb_seq to ensure correct sequence is returned if broken at the first step
                pos_start_pdbseq = pos_start_subseq + res_first + 1
                if pos_start_pdbseq == res_second:
                    return new_pdb_seq
        # Adjust alignment by shifting residues according to connectivity
        gap_matches = re.findall(GAP_ALIGNMENT_PATTERN, subseq)
        if gap_matches:
            residues = "".join(gap_matches)
            residues = residues[rconn_beg:]
            dash_to_add = len(subseq[pos_start_subseq:]) - len(residues)
            new_subseq = subseq[:pos_start_subseq] + "-" * dash_to_add + residues
            new_pdb_seq[res_first + 1 : res_second] = new_subseq
        return new_pdb_seq

    def are_chunks_connected(self, alns):
        """
        Function to checked alignment chunks connectivity and include them in the process if true.

        Args:
            alns (list): list of alignment chunks returned by lalign36. It contains Bio.Align.MultipleSeqAlignment object.

        Returns:
            alignment chunks, only if there are connected and non overlaped.
        """
        # function to check the connectivity of chunks returned by lalign36
        # chunks will be included
        if self.repeated_acc:
            return alns
        alns = [al for al in alns]
        bounds_pdb = [(al[1]._al_start, al[1]._al_stop) for al in alns]
        bounds_unp = [(al[0]._al_start, al[0]._al_stop) for al in alns]
        reference_pdb = bounds_pdb[0]
        reference_unp = bounds_unp[0]
        alns_cleaned = []
        pdb_cleaned = []
        unp_cleaned = []
        for i, bound in enumerate(bounds_pdb):
            if (
                not overlapping([reference_pdb, bound])
                and not overlapping([reference_unp, bounds_unp[i]])
                and not overlapping(
                    pdb_cleaned + [bound]
                )  # check overlap with all included pdb chunks
                and not overlapping(unp_cleaned + [bounds_unp[i]])  # same for unp
            ):
                alns_cleaned.append(alns[i])
                pdb_cleaned.append(bound)
                unp_cleaned.append(bounds_unp[i])
        if alns_cleaned:
            alns_cleaned.insert(0, alns[0])
            connected_chunks = []
            connected_chunks.append(alns_cleaned[0])
            for i in range(len(alns_cleaned) - 1):
                r1 = alns_cleaned[i][1]._al_stop
                r2 = alns_cleaned[i + 1][1]._al_start
                connected = self.check_res_conn(r1, r2)
                if connected:
                    logger.info(
                        f"CONNECTIVITY_CHUNK:{self.chain_obj.pdbid},{self.chain_obj.struct_asym_id},{r1},{r2}"
                    )
                    connected_chunks.append(alns_cleaned[i + 1])
                else:
                    break
            return connected_chunks
        return False

    def alignment_refining(self, al, iso):
        """
        Function to refine extended gaps (i.e more than two dashes) in the alignment from start to end.

        This function aims to move residues to the right continuous region based on connectivity, improving the alignment and continuity between the two regions.

        Args:
            al (Bio.Align.MultipleSeqAlignment object.): the alignment of the chain against one unp accession (iso).
            iso (string): the unp accession.

        Returns:
            al: the new alignment object.
        """
        if self.chain_obj.is_chimera:
            return al
        pdb_seq = str(al[1]._seq)[:]
        start_res = al[1]._al_start
        lenseq = len(pdb_seq)
        if len(pdb_seq) <= 12:
            return al
        new_pdb_seq = list(pdb_seq)
        matches = [
            (match.group(), match.start(), match.end() - 1)
            for match in re.finditer(CONS_ALIGNMENT_PATTERN, pdb_seq)
            if len(match.group()) >= 6
        ]
        if matches:
            handle_beg = True if matches[0][1] != 0 else False
            handle_end = True if matches[-1][-1] != lenseq - 1 else False
            first_cons = matches[0][1]
            last_cons = matches[-1][-1]
            for ind in range(len(matches) - 1):
                res_first = matches[ind][-1]
                res_second = matches[ind + 1][1]
                new_pdb_seq = self.subseq_find_midgap(
                    pdb_seq, new_pdb_seq, res_first, res_second, start_res
                )

            beg_seq = new_pdb_seq[:first_cons]
            last_seq = new_pdb_seq[last_cons + 1 :]
            middle_seq = new_pdb_seq[first_cons : last_cons + 1]
            if handle_beg:
                beg_seq = self.boundaries_check(beg_seq, "start", first_cons, start_res)
            if handle_end:
                last_seq = self.boundaries_check(
                    last_seq, "end", last_cons, al[1]._al_stop
                )
            new_pdb_seq = beg_seq + middle_seq + last_seq

        alignment_seq = "".join(new_pdb_seq)
        if alignment_seq != pdb_seq:
            logger.info(
                f"CONNECTIVITY_REFINING: ENTRY: {self.chain_obj.pdbid}, CHAIN: {self.chain_obj.struct_asym_id},"
                f" ACC: {iso}, PDB_NEW: {alignment_seq}, PDB_OLD: {pdb_seq}, UNP_SEQ: {str(al[0]._seq)[:]}"
            )
        al[1]._seq = Seq(alignment_seq)
        return al

    def check_segments_conn(self):
        """
        Function to check if there are any two segments connected and if true merge them.

        The segments are mapped to the same accession and came from the same chain.

        Returns:
            new_segs: the new list of segments based on connectivity.
        """
        segments = self.chain_obj.segments
        if self.repeated_acc:
            return segments
        new_segs = {k: [] for k in segments}
        for accession in segments:
            segments_acc = segments[accession]
            if len(segments_acc) > 1:
                lst = sorted(segments_acc, key=lambda x: x[0][0])
                current_seg = lst[0]
                for i in range(1, len(lst)):
                    next_seg = lst[i]
                    r1 = current_seg[0][1]
                    r2 = next_seg[0][0]
                    connected = self.check_res_conn(r1, r2)
                    if connected:
                        logger.info(
                            f"CONNECTIVITY_SEG:{self.chain_obj.pdbid},{current_seg},{next_seg}"
                        )
                        new_seg = merge_segment(current_seg, next_seg)
                        current_seg = new_seg
                        new_segs[accession].append(current_seg)
                    else:
                        new_segs[accession].append(current_seg)
                        current_seg = next_seg
                if not connected:
                    new_segs[accession].append(current_seg)
            else:
                return segments
        return new_segs

    def compute_terminals_distance(self):
        """
        Function to compute the distance between two residues and their atom terminals (C and N).

        Returns:
            boolean: True if residues are connected else False.
        """
        if (
            self.r1_ters_coordinates is None
            or self.r2_ters_coordinates is None
            or not self.r1_ters_coordinates["N"]
            or not self.r1_ters_coordinates["C"]
            or not self.r2_ters_coordinates["N"]
            or not self.r2_ters_coordinates["C"]
        ):
            return False
        r1_n = self.r1_ters_coordinates["N"]
        r1_c = self.r1_ters_coordinates["C"]
        r2_n = self.r2_ters_coordinates["N"]
        r2_c = self.r2_ters_coordinates["C"]
        dist1 = compute_atom_distance(r1_n, r2_c)
        dist2 = compute_atom_distance(r1_c, r2_n)
        dist_min = min([dist1, dist2])
        if dist_min <= PEPTIDE_BOND_LEN:
            return True
        return False
