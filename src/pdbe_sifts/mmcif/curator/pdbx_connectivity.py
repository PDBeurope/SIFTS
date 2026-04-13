#!/usr/bin/env python3
"""Connectivity checking and alignment processing for mmCIF structures."""

import regex as re
from gemmi import cif as gemmi_cif

from pdbe_sifts.base.log import logger
from pdbe_sifts.segments_generation.connectivity.ccd_parser import CcdFile
from pdbe_sifts.segments_generation.connectivity.process_connectivity import (
    PEPTIDE_BOND_LEN,
    compute_atom_distance,
)

GAP_ALIGNMENT_PATTERN = r"-+"
CONS_ALIGNMENT_PATTERN = r"[A-Za-z]+"


def apply_fix(gapped_seq: str, p_L: int, p_R: int) -> str:
    """Move the boundary residue at *p_L* to position *p_R* in the gapped string.

    Replaces the residue letter at *p_L* with ``"-"`` and places it at *p_R*
    (the last gap position before the right region, i.e. ``right_start - 1``).

    Example: ``GVT---RVP`` with p_L=2, p_R=5 → ``GV---TRVP``

    Args:
        gapped_seq: Gapped coordinate sequence (letters = observed, ``"-"`` = gap).
        p_L: 0-based index of the boundary residue to move (last of left region).
        p_R: 0-based target position (last ``"-"`` before the right region).

    Returns:
        Corrected gapped string of the same length.
    """
    lst = list(gapped_seq)
    letter = lst[p_L]
    lst[p_L] = "-"
    lst[p_R] = letter
    return "".join(lst)


def apply_multi_fix(
    gapped_seq: str, p_positions: list[int], target_end: int
) -> str:
    """Move multiple boundary residues to just before *target_end*.

    Replaces each source position in *p_positions* with ``"-"`` and writes the
    original letters to the gap positions ending at *target_end*, preserving order.

    Example: ``PPPPPRTPP-----------PPPPP`` with p_positions=[7,8], target_end=19
    → ``PPPPPRT-----------PPPPPPPPP``

    Args:
        gapped_seq: Gapped coordinate sequence string.
        p_positions: Ordered list of source indices to move (left to right).
        target_end: Last target position (= ``p_R - 1``, the last ``"-"`` before
            the right region). The moved letters occupy
            ``[target_end - n + 1 … target_end]``.

    Returns:
        Corrected gapped string of the same length.

    Raises:
        ValueError: If the gap does not contain enough ``"-"`` characters to
            absorb all moved residues.
    """
    n = len(p_positions)
    lst = list(gapped_seq)
    letters = [lst[p] for p in p_positions]
    targets = list(range(target_end - n + 1, target_end + 1))
    if any(lst[t] != "-" for t in targets):
        raise ValueError(
            f"Not enough gap characters to absorb {n} residues at target_end={target_end}"
        )
    for p in p_positions:
        lst[p] = "-"
    for t, letter in zip(targets, letters, strict=False):
        lst[t] = letter
    return "".join(lst)


def get_non_overlapping_chunks(alignments: list) -> list:
    """Return non-overlapping alignment chunks, best-first.

    Keeps the first (highest-score) chunk, then adds each subsequent chunk
    only if it does not overlap — in canonical or coordinate space — with
    any already-kept chunk.  Returns a list containing only the first chunk
    when no non-overlapping secondary chunks are found.

    Args:
        alignments: List of alignment objects exposing
            ``[0]._al_start`` / ``[0]._al_stop`` (canonical positions) and
            ``[1]._al_start`` / ``[1]._al_stop`` (coordinate positions),
            all 1-based inclusive.

    Returns:
        Filtered list of non-overlapping chunks (always at least one element).
        Returns an empty list when *alignments* is empty.
    """
    if not alignments:
        return []
    if len(alignments) == 1:
        return [alignments[0]]

    def _overlaps(a, b):
        return not (a[1] < b[0] or b[1] < a[0])

    kept = [alignments[0]]
    for al in alignments[1:]:
        cb = (al[0]._al_start, al[0]._al_stop)
        qb = (al[1]._al_start, al[1]._al_stop)
        if not any(
            _overlaps(cb, (k[0]._al_start, k[0]._al_stop)) for k in kept
        ) and not any(
            _overlaps(qb, (k[1]._al_start, k[1]._al_stop)) for k in kept
        ):
            kept.append(al)

    return kept


class GapConnectivityChecker:
    """Connectivity checking and alignment processing for a mmCIF structure.

    Designed to handle all polypeptide entity/chain pairs from a single CIF
    file.  Pass the full output of :func:`~pdbe_sifts.seq2seq.run_all` as
    *results* so that every ``(entity_id, chain_id)`` pair is processed.
    """

    def __init__(self, pdbx_mmcif_file, results: list[dict]) -> None:
        """Initialise with a CIF file and seq2seq results.

        Args:
            pdbx_mmcif_file: Path to the mmCIF file.
            alignments: List of result dicts as returned by
                :func:`~pdbe_sifts.seq2seq.run_all` — one dict per
                ``(entity_id, chain_id)`` pair, each containing at minimum
                ``"entity_id"``, ``"chain_id"``, and ``"alignments"`` keys.
        """
        self.cif = pdbx_mmcif_file
        self.results = results
        self.unchunked_alignments: dict[tuple[str, str], list] | None = None

    def process_alignments(self) -> dict[tuple[str, str], list]:
        """Select non-overlapping chunks for every entity/chain pair.

        Iterates over :attr:`results`, calls
        :func:`get_non_overlapping_chunks` for each pair's alignment list,
        and stores the mapping in :attr:`unchunked_alignments`.

        Returns:
            Dict mapping ``(entity_id, chain_id)`` to the list of
            non-overlapping alignment chunks for that pair.  Each list
            contains at least the best (first) chunk.
        """
        self.unchunked_alignments = {}
        for r in self.results:
            key = (r["entity_id"], r["chain_id"])
            self.unchunked_alignments[key] = get_non_overlapping_chunks(
                r["alignments"]
            )
        return self.unchunked_alignments

    def get_gaps_and_cons(self, alignment) -> tuple[list[tuple], list[tuple]]:
        """Extract gap and consecutive regions from a gapped coordinate sequence.

        Args:
            alignment: Single alignment object with ``[1].seq`` (coordinate
                sequence, gapped with ``"-"``).

        Returns:
            ``(gaps, consecutives)`` where each element is a
            ``(substring, start_idx, end_idx)`` tuple with 0-based inclusive
            indices into the gapped string.
        """
        coord_seq = str(alignment[1].seq)
        cons = [
            (match.group(), match.start(), match.end() - 1)
            for match in re.finditer(CONS_ALIGNMENT_PATTERN, coord_seq)
        ]
        gaps = [
            (match.group(), match.start(), match.end() - 1)
            for match in re.finditer(GAP_ALIGNMENT_PATTERN, coord_seq)
        ]
        return gaps, cons

    def _get_ordered_residues(
        self, block, entity_id: str, chain_id: str
    ) -> list[dict]:
        """Build an ordered list of observed residues with all atom coordinates.

        Reads ``_atom_site`` ATOM records for *entity_id* / *chain_id* (model 1
        only) and returns one dict per unique residue, ordered by
        ``(auth_seq_id, ins_code)``.  Each dict stores every atom found for
        that residue, keyed by ``(label_atom_id, ordinal)`` to match the
        format expected by :class:`~pdbe_sifts.segments_generation.connectivity.ccd_parser.CcdFile`.

        Args:
            block: :class:`gemmi.cif.Block` of the mmCIF file.
            entity_id: Entity identifier string (e.g. ``"1"``).
            chain_id: Author asymmetric-unit chain identifier (e.g. ``"A"``).

        Returns:
            List of dicts ordered by residue position, each with keys:

            * ``auth_seq_id`` (:class:`int`) — author sequence number.
            * ``ins_code`` (:class:`str`) — insertion code (``"."`` when absent).
            * ``comp_id`` (:class:`str`) — three-letter component identifier.
            * ``atoms`` (:class:`dict`) — maps ``(label_atom_id, ordinal_int)``
              to ``[x, y, z]`` float coordinate lists.
        """
        atom_site = block.get_mmcif_category("_atom_site")
        if not atom_site:
            return []

        groups = atom_site.get("group_PDB", [])
        entity_ids = atom_site.get("label_entity_id", [])
        chain_ids = atom_site.get("auth_asym_id", [])
        seq_ids = atom_site.get("auth_seq_id", [])
        ins_codes = atom_site.get("pdbx_PDB_ins_code", [])
        comp_ids = atom_site.get("label_comp_id", [])
        atom_ids = atom_site.get("label_atom_id", [])
        xs = atom_site.get("Cartn_x", [])
        ys = atom_site.get("Cartn_y", [])
        zs = atom_site.get("Cartn_z", [])
        models = atom_site.get("pdbx_PDB_model_num", [])

        # Collect atoms per residue key, tracking ordinals per residue
        residue_order: list[
            tuple
        ] = []  # ordered unique (auth_seq_id, ins_code) keys
        residue_data: dict[tuple, dict] = {}
        prev_key: tuple | None = None
        atom_ordinal = 1

        for i, group in enumerate(groups):
            # Accept HETATM too: modified residues (e.g. HYP, MSE) belonging
            # to a polypeptide entity are often tagged HETATM in the CIF.
            if group not in ("ATOM", "HETATM"):
                continue
            if entity_ids[i] != entity_id:
                continue
            if chain_ids[i] != chain_id:
                continue
            if models and models[i] not in ("1", 1):
                continue

            ins = ins_codes[i] if ins_codes else "."
            key = (int(seq_ids[i]), ins)

            if key != prev_key:
                atom_ordinal = 1
                prev_key = key
                if key not in residue_data:
                    residue_order.append(key)
                    residue_data[key] = {
                        "auth_seq_id": int(seq_ids[i]),
                        "ins_code": ins,
                        "comp_id": comp_ids[i],
                        "atoms": {},
                    }

            atom_key = (atom_ids[i], atom_ordinal)
            residue_data[key]["atoms"][atom_key] = [
                float(xs[i]),
                float(ys[i]),
                float(zs[i]),
            ]
            atom_ordinal += 1

        return [residue_data[k] for k in residue_order]

    def _get_terminal_coords(
        self, residue: dict, terminal: str
    ) -> list[float] | None:
        """Return C or N backbone atom coordinates for a residue using CcdFile.

        Looks up the terminal atom label and ordinal from the Chemical Component
        Dictionary, then retrieves the corresponding coordinates from the
        residue's ``"atoms"`` dict.

        Args:
            residue: Residue dict as returned by :meth:`_get_ordered_residues`.
            terminal: ``"C"`` for the C-terminal backbone atom or ``"N"`` for
                the N-terminal backbone atom.

        Returns:
            ``[x, y, z]`` float list, or ``None`` if the CCD lookup fails or
            the atom is absent from the structure.
        """
        ters = CcdFile(residue["comp_id"]).process()
        if ters == 0:
            logger.warning(
                f"CCD lookup failed for {residue['comp_id']}, cannot get {terminal} coords"
            )
            return None
        atom_label = ters[terminal][0]
        atom_ordinal = int(ters[terminal][1])
        return residue["atoms"].get((atom_label, atom_ordinal))

    def _check_continuity(
        self,
        ordered: list[dict],
        start_idx: int,
        direction: str,
        steps: int = 10,
    ) -> dict:
        """Walk the residue chain checking C→N distances at each step.

        Args:
            ordered: Residue list from :meth:`_get_ordered_residues`.
            start_idx: 0-based index to start from.
            direction: ``"left"`` to walk backwards, ``"right"`` to walk forwards.
            steps: Maximum number of consecutive bonds to check.

        Returns:
            Dict with ``"distances"`` (list of float or None per step) and
            ``"first_break"`` (0-based step index of the first break, or None
            if the chain is fully connected within *steps*).
        """
        distances = []
        first_break = None
        for i in range(steps):
            if direction == "left":
                curr, prev = start_idx - i, start_idx - i - 1
                if prev < 0 or curr >= len(ordered):
                    break
                c = self._get_terminal_coords(ordered[prev], "C")
                n = self._get_terminal_coords(ordered[curr], "N")
            else:
                curr, nxt = start_idx + i, start_idx + i + 1
                if nxt >= len(ordered):
                    break
                c = self._get_terminal_coords(ordered[curr], "C")
                n = self._get_terminal_coords(ordered[nxt], "N")
            if c is None or n is None:
                distances.append(None)
                if first_break is None:
                    first_break = i
            else:
                d = compute_atom_distance(c, n)
                distances.append(d)
                if d > PEPTIDE_BOND_LEN and first_break is None:
                    first_break = i
        return {"distances": distances, "first_break": first_break}

    def _check_single_gap(
        self,
        gapped: str,
        ordered: list[dict],
        base: int,
        gap_bounds: tuple[str, int, int],
        cons: list[tuple],
    ) -> dict | None:
        """Analyse a single gap and return the info dict, or None for border gaps.

        Computes the two key distances around the boundary residue ``rL`` (last of
        the left observed region) and decides whether it belongs to the left block,
        the right block, or is ambiguous.  When a fix is warranted, stores both the
        pre-computed ``fixed_seq`` (for display) and the raw ``source_positions`` /
        ``target_end`` parameters so that callers can re-apply the same fix on an
        evolving string without overwriting earlier fixes.

        Args:
            gapped: Current gapped coordinate string (letters = observed, ``"-"`` = gap).
            ordered: Residue list from :meth:`_get_ordered_residues` (unchanged across
                iterations for the same chunk).
            base: ``chunk[1]._al_start - 1`` — number of coordinate residues that
                precede this chunk; used to convert gapped positions to ``ordered``
                indices.
            gap_bounds: ``(substring, gap_start, gap_end)`` tuple (0-based inclusive)
                for the gap being analysed.
            cons: List of ``(substring, start, end)`` tuples for all consecutive
                (non-gap) regions in *gapped*, as returned by the regex scanner.

        Returns:
            Info dict with keys ``gap_pos``, ``rL``, ``rL_num``, ``rL_prev``,
            ``rL_prev_num``, ``rR``, ``rR_num``, ``d_conn_L``, ``d_conn_R``,
            ``action``, ``original_seq``, and optionally ``fixed_seq``,
            ``source_positions``, ``target_end``, ``left_continuity``,
            ``right_continuity``.  Returns ``None`` when the gap is at the edge
            of the alignment (no left or right region) or when a residue index
            is out of range.
        """
        _, gap_start, gap_end = gap_bounds

        left_r = next((c for c in cons if c[2] == gap_start - 1), None)
        right_r = next((c for c in cons if c[1] == gap_end + 1), None)
        if left_r is None or right_r is None:
            return None

        p_L = left_r[2]
        p_R = right_r[1]
        p_L_prev = p_L - 1 if p_L - 1 >= left_r[1] else None

        def _idx(pos: int) -> int:
            return base + sum(1 for c in gapped[: pos + 1] if c != "-") - 1

        idx_rL = _idx(p_L)
        idx_rR = _idx(p_R)
        idx_rL_prev = _idx(p_L_prev) if p_L_prev is not None else None

        if idx_rL >= len(ordered) or idx_rR >= len(ordered):
            logger.warning(f"Residue index out of range at gap ({p_L}, {p_R})")
            return None

        res_rL = ordered[idx_rL]
        res_rR = ordered[idx_rR]
        res_rL_prev = (
            ordered[idx_rL_prev]
            if idx_rL_prev is not None and idx_rL_prev >= 0
            else None
        )

        c_prev = (
            self._get_terminal_coords(res_rL_prev, "C") if res_rL_prev else None
        )
        n_rL = self._get_terminal_coords(res_rL, "N")
        c_rL = self._get_terminal_coords(res_rL, "C")
        n_rR = self._get_terminal_coords(res_rR, "N")

        d_conn_L = (
            compute_atom_distance(c_prev, n_rL) if c_prev and n_rL else None
        )
        d_conn_R = compute_atom_distance(c_rL, n_rR) if c_rL and n_rR else None

        conn_L = d_conn_L is not None and d_conn_L <= PEPTIDE_BOND_LEN
        conn_R = d_conn_R is not None and d_conn_R <= PEPTIDE_BOND_LEN

        info: dict = {
            "gap_pos": (p_L, p_R),
            "rL": gapped[p_L],
            "rL_num": res_rL["auth_seq_id"],
            "rL_prev": gapped[p_L_prev] if p_L_prev is not None else None,
            "rL_prev_num": res_rL_prev["auth_seq_id"] if res_rL_prev else None,
            "rR": gapped[p_R],
            "rR_num": res_rR["auth_seq_id"],
            "d_conn_L": d_conn_L,
            "d_conn_R": d_conn_R,
            "original_seq": gapped,
        }

        if conn_L and not conn_R:
            info["action"] = "ok"

        elif conn_R and not conn_L:
            info["action"] = "fix"
            info["source_positions"] = [p_L]
            info["target_end"] = p_R - 1
            info["fixed_seq"] = apply_fix(gapped, p_L, p_R - 1)

        elif conn_L and conn_R:
            left_start = idx_rL_prev if idx_rL_prev is not None else idx_rL
            left_cont = self._check_continuity(ordered, left_start, "left")
            right_cont = self._check_continuity(ordered, idx_rR, "right")
            left_break = left_cont["first_break"]
            right_break = right_cont["first_break"]
            info["left_continuity"] = left_cont
            info["right_continuity"] = right_cont

            if right_break is not None and left_break is None:
                info["action"] = "ok"
            elif left_break is not None and right_break is None:
                if left_break == 0 and p_L_prev is not None:
                    info["action"] = "fix"
                    info["source_positions"] = [p_L_prev, p_L]
                    info["target_end"] = p_R - 1
                    info["fixed_seq"] = apply_multi_fix(
                        gapped, [p_L_prev, p_L], p_R - 1
                    )
                else:
                    info["action"] = "fix_partial"
                    info["source_positions"] = [p_L]
                    info["target_end"] = p_R - 1
                    info["fixed_seq"] = apply_fix(gapped, p_L, p_R - 1)
            else:
                info["action"] = "ambiguous"

        else:
            info["action"] = "chain_break"

        return info

    def check_gap_connectivity(self) -> dict[tuple[str, str], list[dict]]:
        """Check C→N peptide bond connectivity at every inter-region gap boundary.

        For each gap separating two consecutive observed regions, tests whether
        the boundary residue ``rL`` (last of the left region) belongs to the
        left block or the right block by measuring two distances:

        * ``d_conn_L`` : C(rL_prev) → N(rL) — rL follows rL_prev in the chain
        * ``d_conn_R`` : C(rL) → N(rR) — rL precedes rR in the chain

        Decision table (threshold :data:`PEPTIDE_BOND_LEN` = 1.42 Å):

        =========  =========  ============================================================
        conn_L     conn_R     action
        =========  =========  ============================================================
        True       False      ``"ok"`` — rL belongs to the left block
        False      True       ``"fix"`` — rL moves just before rR
        True       True       extended check (±10 residues); ``"ok"``, ``"fix"``,
                              or ``"ambiguous"``
        False      False      ``"chain_break"``
        =========  =========  ============================================================

        Requires :meth:`process_alignments` to have been called first.

        Returns:
            Dict mapping ``(entity_id, chain_id)`` to a list of per-gap dicts.
            Each dict contains: ``chunk_idx``, ``gap_pos``, ``rL``, ``rL_prev``,
            ``rR``, ``d_conn_L``, ``d_conn_R``, ``action``, ``original_seq``,
            and optionally ``fixed_seq`` (when action is ``"fix"``),
            ``source_positions`` and ``target_end`` (parameters to re-apply the
            fix on any string), ``left_continuity`` and ``right_continuity``
            (when both connections are detected).

        Raises:
            RuntimeError: If :meth:`process_alignments` has not been called.
        """
        if self.unchunked_alignments is None:
            raise RuntimeError(
                "Call process_alignments() before check_gap_connectivity()."
            )

        block = gemmi_cif.read(str(self.cif)).sole_block()
        output: dict[tuple[str, str], list[dict]] = {}

        for r in self.results:
            entity_id, chain_id = r["entity_id"], r["chain_id"]
            key = (entity_id, chain_id)
            ordered = self._get_ordered_residues(block, entity_id, chain_id)
            infos: list[dict] = []

            for chunk_idx, chunk in enumerate(
                self.unchunked_alignments.get(key, [])
            ):
                gapped = str(chunk[1].seq)
                gaps, cons = self.get_gaps_and_cons(chunk)
                base = chunk[1]._al_start - 1

                for gap_bounds in gaps:
                    info = self._check_single_gap(
                        gapped, ordered, base, gap_bounds, cons
                    )
                    if info is None:
                        continue
                    info["chunk_idx"] = chunk_idx

                    # Log context-rich messages now that we have entity/chain info
                    action = info["action"]
                    p_L, p_R = info["gap_pos"]
                    if action == "ok" and "left_continuity" in info:
                        rb = info["right_continuity"]["first_break"]
                        logger.info(
                            f"{entity_id}/{chain_id}: rL={info['rL']} belongs to left "
                            f"block (right chain break at step {rb})"
                        )
                    elif action == "fix" and "left_continuity" in info:
                        lb = info["left_continuity"]["first_break"]
                        logger.info(
                            f"{entity_id}/{chain_id}: rL_prev+rL moving to right block "
                            f"(left break at step {lb})"
                        )
                    elif action == "fix_partial":
                        lb = info.get("left_continuity", {}).get("first_break")
                        logger.warning(
                            f"{entity_id}/{chain_id}: left break at step {lb}, "
                            "partial fix applied — manual review may be needed"
                        )
                    elif action == "ambiguous":
                        logger.warning(
                            f"{entity_id}/{chain_id}: ambiguous connectivity at "
                            f"gap ({p_L}, {p_R})"
                        )

                    infos.append(info)

            output[key] = infos

        return output

    def show(self, context: int = 8) -> None:
        """Print a human-readable summary of gap connectivity results.

        Calls :meth:`check_gap_connectivity` if not already cached, then
        prints one block per ``(entity_id, chain_id)`` pair with, for each
        gap: the action, the boundary residues with their author sequence
        numbers, the C→N distances, and (when a fix is applied) the alignment
        string before and after correction around the gap.

        Args:
            context: Number of characters to show on each side of the gap in
                the before/after alignment snippet.
        """
        gaps = self.check_gap_connectivity()

        action_label = {
            "ok": "OK          ",
            "fix": "FIX         ",
            "fix_partial": "FIX PARTIAL ",
            "chain_break": "CHAIN BREAK ",
            "ambiguous": "AMBIGUOUS   ",
        }

        for (entity_id, chain_id), infos in gaps.items():
            print(
                f"\nentity={entity_id}  chain={chain_id}  ({len(infos)} gap(s))"
            )
            if not infos:
                print("  (no inter-region gaps)")
                continue

            for info in infos:
                action = info["action"]
                label = action_label.get(action, f"{action:12}")

                prev_str = (
                    f"{info['rL_prev']}({info['rL_prev_num']}) - "
                    if info["rL_prev"] is not None
                    else ""
                )
                d_L = (
                    f"{info['d_conn_L']:.2f}"
                    if info["d_conn_L"] is not None
                    else "n/a"
                )
                d_R = (
                    f"{info['d_conn_R']:.2f}"
                    if info["d_conn_R"] is not None
                    else "n/a"
                )

                print(
                    f"  [{label}]  "
                    f"{prev_str}{info['rL']}({info['rL_num']}) ~~gap~~ "
                    f"{info['rR']}({info['rR_num']})   "
                    f"d_L={d_L} Å  d_R={d_R} Å"
                )

                if action in ("fix", "fix_partial"):
                    orig = info["original_seq"]
                    fixed = info["fixed_seq"]
                    p_L, p_R = info["gap_pos"]
                    lo = max(0, p_L - context)
                    hi = min(len(orig), p_R + context + 1)
                    print(f"    before : ...{orig[lo:hi]}...")
                    print(f"    after  : ...{fixed[lo:hi]}...")

                if (
                    action in ("ambiguous", "fix", "fix_partial")
                    and "left_continuity" in info
                ):
                    lc = info["left_continuity"]
                    rc = info["right_continuity"]
                    print(
                        f"    left  first_break={lc['first_break']}  distances={lc['distances']}"
                    )
                    print(
                        f"    right first_break={rc['first_break']}  distances={rc['distances']}"
                    )

    def get_corrected_sequences(
        self, max_iter: int = 10
    ) -> dict[tuple[str, str], list[dict]]:
        """Apply connectivity fixes iteratively until each chunk's sequence is stable.

        For each alignment chunk, repeatedly re-analyses all inter-region gaps on
        the current (possibly already corrected) gapped string, applies every
        ``"fix"`` / ``"fix_partial"`` action found, and stops when no further fixes
        are produced or *max_iter* passes have been completed.

        Because :func:`apply_fix` and :func:`apply_multi_fix` are pure character
        swaps (string length is constant, no positions ever shift), multiple fixes
        within a single pass operate on non-overlapping segments and may be applied
        in any order with identical results.  The cascade effect — where fixing one
        gap changes the boundary context of a neighbouring gap — is handled by
        re-running the full gap analysis at the start of each pass.

        Requires :meth:`process_alignments` to have been called first.

        Args:
            max_iter: Maximum number of correction passes per chunk before
                giving up.  Defaults to ``10``.

        Returns:
            Dict mapping ``(entity_id, chain_id)`` to a list of per-chunk dicts,
            each containing:

            * ``chunk_idx`` (:class:`int`) — index of the alignment chunk.
            * ``original``  (:class:`str`) — gapped coordinate string before any fix.
            * ``corrected`` (:class:`str`) — final string after all cascaded fixes.
            * ``n_iter``    (:class:`int`) — number of passes that produced at
              least one fix (``0`` means no fixes were needed).

        Raises:
            RuntimeError: If :meth:`process_alignments` has not been called.
        """
        if self.unchunked_alignments is None:
            raise RuntimeError(
                "Call process_alignments() before get_corrected_sequences()."
            )

        block = gemmi_cif.read(str(self.cif)).sole_block()
        output: dict[tuple[str, str], list[dict]] = {}

        for r in self.results:
            entity_id, chain_id = r["entity_id"], r["chain_id"]
            key = (entity_id, chain_id)
            ordered = self._get_ordered_residues(block, entity_id, chain_id)
            chunk_results: list[dict] = []

            for chunk_idx, chunk in enumerate(
                self.unchunked_alignments.get(key, [])
            ):
                original = str(chunk[1].seq)
                current = original
                base = chunk[1]._al_start - 1
                n_iter = 0

                for _ in range(max_iter):
                    # Re-detect gaps and consecutive regions on the current string
                    cons = [
                        (m.group(), m.start(), m.end() - 1)
                        for m in re.finditer(CONS_ALIGNMENT_PATTERN, current)
                    ]
                    gaps = [
                        (m.group(), m.start(), m.end() - 1)
                        for m in re.finditer(GAP_ALIGNMENT_PATTERN, current)
                    ]

                    fixes = []
                    for gap_bounds in gaps:
                        info = self._check_single_gap(
                            current, ordered, base, gap_bounds, cons
                        )
                        if info is not None and info["action"] in (
                            "fix",
                            "fix_partial",
                        ):
                            fixes.append(info)

                    if not fixes:
                        break  # sequence is stable

                    # Apply each fix via stored parameters rather than pre-computed
                    # fixed_seq: all fixed_seq values were computed on the same
                    # 'current', so assigning them sequentially would make each
                    # one overwrite the previous.  Using source_positions +
                    # target_end re-applies the same logical move on the evolving
                    # string instead.
                    for fix in fixes:
                        src = fix["source_positions"]
                        tgt = fix["target_end"]
                        if len(src) == 1:
                            current = apply_fix(current, src[0], tgt)
                        else:
                            current = apply_multi_fix(current, src, tgt)
                    n_iter += 1

                chunk_results.append(
                    {
                        "chunk_idx": chunk_idx,
                        "original": original,
                        "corrected": current,
                        "n_iter": n_iter,
                    }
                )

            output[key] = chunk_results

        return output

    def get_merged_alignment(self) -> dict[tuple[str, str], dict]:
        """Return one merged alignment per (entity_id, chain_id).

        Orchestrates the full connectivity pipeline:

        1. :meth:`process_alignments` — select non-overlapping chunks.
        2. :meth:`get_corrected_sequences` — correct gap boundaries.
        3. Merge all corrected chunks into a single canonical-length pair.

        Each chunk's corrected coordinate sequence is projected back onto the
        canonical sequence using the chunk's ``_al_start`` position and the
        canonical gapped sequence to track alignment columns.  Canonical
        positions not covered by any chunk remain ``"-"``.

        Returns:
            Dict mapping ``(entity_id, chain_id)`` to:

            * ``"canonical"`` (:class:`str`) — full canonical sequence, no gaps.
            * ``"coord"``     (:class:`str`) — same length; ``"-"`` where unobserved.
            * ``"n_fixes"``   (:class:`int`) — chunks where corrections were applied.
        """
        self.process_alignments()
        corrected = self.get_corrected_sequences()
        merged: dict[tuple[str, str], dict] = {}

        for r in self.results:
            key = (r["entity_id"], r["chain_id"])
            canonical = r["canonical"]
            chunks = self.unchunked_alignments.get(key, [])
            chunk_corrections = {
                c["chunk_idx"]: c["corrected"] for c in corrected.get(key, [])
            }

            coord_merged = ["-"] * len(canonical)
            for i, chunk in enumerate(chunks):
                coord_gapped = chunk_corrections.get(i, str(chunk[1].seq))
                can_gapped = str(chunk[0].seq)
                can_pos = chunk[0]._al_start - 1  # 0-based index into canonical

                for can_ch, coord_ch in zip(
                    can_gapped, coord_gapped, strict=False
                ):
                    if can_ch != "-":
                        if coord_ch != "-":
                            coord_merged[can_pos] = coord_ch
                        can_pos += 1

            n_fixes = sum(
                1
                for c in corrected.get(key, [])
                if c["corrected"] != c["original"]
            )
            merged[key] = {
                "canonical": canonical,
                "coord": "".join(coord_merged),
                "n_fixes": n_fixes,
            }

        return merged
