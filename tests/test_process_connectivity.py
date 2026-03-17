"""Tests for pdbe_sifts.segments_generation.connectivity.process_connectivity."""

from types import SimpleNamespace
from unittest.mock import patch

import pytest
from pdbe_sifts.segments_generation.connectivity.process_connectivity import (
    ConnectivityCheck,
    compute_atom_distance,
    fmt_ranges,
    is_valid_residue_status,
    merge_segment,
    overlapping,
)

PEPTIDE_BOND_LEN = 1.42


# ── Pure function tests ───────────────────────────────────────────────────────


class TestFmtRanges:
    def test_string_ranges(self):
        assert fmt_ranges(["1-10", "20-30"]) == [[1, 10], [20, 30]]

    def test_list_ranges_unchanged(self):
        assert fmt_ranges([[1, 10], [20, 30]]) == [[1, 10], [20, 30]]

    def test_mixed(self):
        assert fmt_ranges(["5-15", [20, 25]]) == [[5, 15], [20, 25]]


class TestOverlapping:
    def test_non_overlapping(self):
        assert overlapping([[1, 10], [11, 20]]) is False

    def test_overlapping(self):
        assert overlapping([[1, 15], [10, 20]]) is True

    def test_adjacent_not_overlapping(self):
        # 10 is end of first, 10 is start of second → overlapping by ≤ check
        assert overlapping([[1, 10], [10, 20]]) is True

    def test_single_range(self):
        assert overlapping([[1, 100]]) is False

    def test_contained(self):
        assert overlapping([[1, 100], [10, 50]]) is True


class TestMergeSegment:
    def test_basic_merge(self):
        seg1 = ((1, 50), (1, 50))
        seg2 = ((51, 100), (51, 100))
        result = merge_segment(seg1, seg2)
        assert result == ((1, 100), (1, 100))

    def test_pdb_bounds(self):
        seg1 = ((10, 30), (1, 20))
        seg2 = ((35, 60), (25, 50))
        pdb_start, pdb_end = merge_segment(seg1, seg2)[0]
        assert pdb_start == 10
        assert pdb_end == 60

    def test_unp_bounds(self):
        seg1 = ((10, 30), (1, 20))
        seg2 = ((35, 60), (25, 50))
        unp_start, unp_end = merge_segment(seg1, seg2)[1]
        assert unp_start == 1
        assert unp_end == 50


class TestComputeAtomDistance:
    def test_zero_distance(self):
        assert compute_atom_distance([0, 0, 0], [0, 0, 0]) == 0.0

    def test_unit_distance_x(self):
        d = compute_atom_distance([0, 0, 0], [1, 0, 0])
        assert d == 1.0

    def test_3d_distance(self):
        # sqrt(3) ≈ 1.73
        d = compute_atom_distance([0, 0, 0], [1, 1, 1])
        assert 1.73 <= d <= 1.74

    def test_truncation(self):
        # 1.41421... → truncated to 1.41
        d = compute_atom_distance([0, 0, 0], [1, 1, 0])
        assert d == 1.41

    def test_symmetry(self):
        a, b = [1, 2, 3], [4, 5, 6]
        assert compute_atom_distance(a, b) == compute_atom_distance(b, a)


class TestIsValidResidueStatus:
    def _r(self, rtype="ATOM", observed=True, auth_ins=False):
        return SimpleNamespace(
            rtype=rtype, observed=observed, auth_ins=auth_ins
        )

    def test_valid_pair(self):
        assert is_valid_residue_status(self._r(), self._r()) is True

    @pytest.mark.parametrize("rtype", ["Insertion", "Linker", "Chromophore"])
    def test_forbidden_rtype_r1(self, rtype):
        assert is_valid_residue_status(self._r(rtype=rtype), self._r()) is False

    @pytest.mark.parametrize("rtype", ["Insertion", "Linker", "Chromophore"])
    def test_forbidden_rtype_r2(self, rtype):
        assert is_valid_residue_status(self._r(), self._r(rtype=rtype)) is False

    def test_unobserved_r1(self):
        assert (
            is_valid_residue_status(self._r(observed=False), self._r()) is False
        )

    def test_unobserved_r2(self):
        assert (
            is_valid_residue_status(self._r(), self._r(observed=False)) is False
        )

    def test_insertion_code_r1(self):
        assert (
            is_valid_residue_status(self._r(auth_ins="A"), self._r()) is False
        )

    def test_insertion_code_r2(self):
        assert (
            is_valid_residue_status(self._r(), self._r(auth_ins="B")) is False
        )


# ── ConnectivityCheck.check_segments_conn tests ───────────────────────────────


def _make_chain(segments, residues=None):
    return SimpleNamespace(
        segments=segments,
        residues=residues or [],
        pdbid="1TST",
        is_chimera=False,
    )


class TestCheckSegmentsConn:
    """Mirror the tests in test_connectivity_claude but for the original module."""

    def _check(self, segments, connected_pairs=None):
        """Run check_segments_conn with mocked check_res_conn."""
        connected_pairs = connected_pairs or set()
        chain = _make_chain(segments)
        cc = ConnectivityCheck(chain_obj=chain)

        def mock_check(r1, r2):
            return (r1, r2) in connected_pairs

        with patch.object(cc, "check_res_conn", side_effect=mock_check):
            return cc.check_segments_conn()

    def test_single_segment_returned_unchanged(self):
        segs = {"P00001": [((1, 50), (1, 50))]}
        result = self._check(segs)
        assert result == segs

    def test_two_disconnected_segments(self):
        segs = {"P00001": [((1, 30), (1, 30)), ((50, 80), (50, 80))]}
        result = self._check(segs, connected_pairs=set())
        assert len(result["P00001"]) == 2

    def test_two_connected_segments_merge(self):
        segs = {"P00001": [((1, 30), (1, 30)), ((31, 80), (31, 80))]}
        result = self._check(segs, connected_pairs={(30, 31)})
        assert len(result["P00001"]) == 1
        assert result["P00001"][0] == ((1, 80), (1, 80))

    def test_three_all_connected(self):
        segs = {
            "P00001": [
                ((1, 20), (1, 20)),
                ((21, 40), (21, 40)),
                ((41, 60), (41, 60)),
            ]
        }
        result = self._check(segs, connected_pairs={(20, 21), (40, 41)})
        assert len(result["P00001"]) == 1
        assert result["P00001"][0] == ((1, 60), (1, 60))

    def test_three_first_two_connected_last_not(self):
        segs = {
            "P00001": [
                ((1, 20), (1, 20)),
                ((21, 40), (21, 40)),
                ((50, 70), (50, 70)),
            ]
        }
        result = self._check(segs, connected_pairs={(20, 21)})
        assert len(result["P00001"]) == 2
        merged, last = result["P00001"]
        assert merged == ((1, 40), (1, 40))
        assert last == ((50, 70), (50, 70))

    def test_no_double_append_when_connected(self):
        """Three segments all connected must yield exactly 1 merged segment."""
        segs = {
            "P00001": [
                ((1, 10), (1, 10)),
                ((11, 20), (11, 20)),
                ((21, 30), (21, 30)),
            ]
        }
        result = self._check(segs, connected_pairs={(10, 11), (20, 21)})
        assert len(result["P00001"]) == 1

    def test_repeated_acc_returns_segments_unchanged(self):
        segs = {"P00001": [((1, 30), (1, 30)), ((31, 80), (31, 80))]}
        chain = _make_chain(segs)
        cc = ConnectivityCheck(chain_obj=chain, repeated_acc=True)
        result = cc.check_segments_conn()
        assert result == segs

    def test_multiple_accessions_independent(self):
        segs = {
            "P00001": [((1, 30), (1, 30)), ((31, 60), (31, 60))],
            "P00002": [((1, 50), (1, 50)), ((60, 100), (60, 100))],
        }
        result = self._check(segs, connected_pairs={(30, 31)})
        assert len(result["P00001"]) == 1
        assert len(result["P00002"]) == 2
