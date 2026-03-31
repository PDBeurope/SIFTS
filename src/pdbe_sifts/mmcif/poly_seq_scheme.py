"""Reconstruct ``_pdbx_poly_seq_scheme`` from ``_entity_poly_seq`` and ``_atom_site``.

Some CIF files (in-house, experimental, or partial exports) lack the
``_pdbx_poly_seq_scheme`` category that the SIFTS pipeline requires.  This
module provides a self-contained reconstruction pipeline:

1. Read the canonical sequence from ``_entity_poly_seq`` (three-letter codes,
   one row per position).  Falls back to ``_entity_poly`` one-letter fields
   (``pdbx_seq_one_letter_code_can``, then ``pdbx_seq_one_letter_code``)
   when ``_entity_poly_seq`` is absent.
2. Read the observed residues from ``_atom_site`` ATOM records (model 1 only).
3. Align the two one-letter sequences with ``lalign36`` to establish the
   correspondence between canonical positions and coordinate residues.
4. Return a dict of parallel lists that can be written back to a gemmi block
   with ``block.set_mmcif_category("_pdbx_poly_seq_scheme", scheme)``.
"""

import logging

from pdbe_sifts.mmcif.chem_comp import STANDARD_AA, ChemCompMapping
from pdbe_sifts.segments_generation.alignment import do_alignment_lalign36

logger = logging.getLogger(__name__)

_MISSING = "?"

# Reverse mapping: one-letter → three-letter (standard amino acids only)
_ONE_TO_THREE: dict[str, str] = {v: k for k, v in STANDARD_AA.items()}

# Ordered column names for _pdbx_poly_seq_scheme
SCHEME_COLUMNS: tuple[str, ...] = (
    "asym_id",
    "entity_id",
    "seq_id",
    "mon_id",
    "ndb_seq_num",
    "pdb_seq_num",
    "auth_seq_num",
    "pdb_mon_id",
    "auth_mon_id",
    "pdb_strand_id",
    "pdb_ins_code",
    "hetero",
)


def get_canonical_residues(block, entity_id: str) -> list[tuple[int, str]]:
    """Return the canonical residue list for *entity_id*, with fallback.

    Tries ``_entity_poly_seq`` first (preferred — three-letter codes per
    position).  When that category is absent or yields nothing for
    *entity_id*, falls back to ``_entity_poly`` one-letter fields in priority
    order: ``pdbx_seq_one_letter_code_can`` (modified residues normalised to
    their standard equivalents) then ``pdbx_seq_one_letter_code`` (raw form).
    Each one-letter code is back-converted to a standard three-letter code via
    :data:`_ONE_TO_THREE`; unknown codes map to ``"UNK"``.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        entity_id: Entity identifier string (e.g. ``"1"``).

    Returns:
        List of ``(seq_id, three_letter_code)`` tuples sorted by *seq_id*.
        Returns an empty list when neither source contains the entity.
    """
    # Primary: _entity_poly_seq
    eps = block.get_mmcif_category("_entity_poly_seq")
    if eps:
        result: dict[int, str] = {}
        for idx, eid in enumerate(eps.get("entity_id", [])):
            if eid != entity_id:
                continue
            num = int(eps["num"][idx])
            if num not in result:
                result[num] = eps["mon_id"][idx]
        if result:
            return sorted(result.items())

    # Fallback: _entity_poly.pdbx_seq_one_letter_code_can, then pdbx_seq_one_letter_code
    logger.debug(
        f"_entity_poly_seq not available for entity {entity_id}, "
        "falling back to _entity_poly one-letter sequence fields"
    )
    entity_poly = block.get_mmcif_category("_entity_poly")
    if entity_poly:
        for idx, eid in enumerate(entity_poly.get("entity_id", [])):
            if eid != entity_id:
                continue
            # Try canonical form first, then raw one-letter code
            for field in (
                "pdbx_seq_one_letter_code_can",
                "pdbx_seq_one_letter_code",
            ):
                seq_list = entity_poly.get(field, [])
                if idx >= len(seq_list):
                    continue
                raw = seq_list[idx]
                if raw:
                    seq_1L = raw.replace("\n", "").replace(" ", "")
                    return [
                        (i + 1, _ONE_TO_THREE.get(aa, "UNK"))
                        for i, aa in enumerate(seq_1L)
                    ]

    return []


def get_coordinate_residues(
    block, entity_id: str, chain_id: str
) -> tuple[list[dict], str]:
    """Return observed residues for *entity_id* / *chain_id* from ``_atom_site``.

    Reads ATOM records from model 1, deduplicates by ``(auth_seq_id,
    pdbx_PDB_ins_code)``, and detects microheterogeneity.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        entity_id: Entity identifier string (e.g. ``"1"``).
        chain_id: Author chain identifier (e.g. ``"A"``).

    Returns:
        A ``(residues, asym_id)`` tuple where:

        * *residues* is a list of dicts — one per unique residue — ordered by
          ``(auth_seq_id, ins_code)``, each containing keys ``auth_seq_id``,
          ``ins_code``, ``comp_id``, ``asym_id``, and ``hetero``.
        * *asym_id* is the ``label_asym_id`` shared by all residues in the
          chain (used to fill unobserved rows); ``""`` when nothing is found.
    """
    atom_site = block.get_mmcif_category("_atom_site")
    if not atom_site:
        return [], ""

    groups = atom_site.get("group_PDB", [])
    entity_ids = atom_site.get("label_entity_id", [])
    chain_ids = atom_site.get("auth_asym_id", [])
    seq_ids = atom_site.get("auth_seq_id", [])
    ins_codes = atom_site.get("pdbx_PDB_ins_code", [])
    comp_ids = atom_site.get("label_comp_id", [])
    asym_ids = atom_site.get("label_asym_id", [])
    model_nums = atom_site.get("pdbx_PDB_model_num", [])

    seen_comp: dict[tuple, str] = {}
    hetero_pos: set[tuple] = set()
    seen_asym: str = ""

    for i, group in enumerate(groups):
        if group != "ATOM":
            continue
        if entity_ids[i] != entity_id:
            continue
        if chain_ids[i] != chain_id:
            continue
        if model_nums and model_nums[i] not in ("1", 1):
            continue

        ins = ins_codes[i] if ins_codes else "."
        key = (int(seq_ids[i]), ins)
        comp = comp_ids[i]

        if key not in seen_comp:
            seen_comp[key] = comp
            if asym_ids:
                seen_asym = asym_ids[i]
        elif seen_comp[key] != comp:
            hetero_pos.add(key)

    residues = [
        {
            "auth_seq_id": k[0],
            "ins_code": k[1],
            "comp_id": v,
            "asym_id": seen_asym,
            "hetero": "y" if k in hetero_pos else "n",
        }
        for k, v in sorted(seen_comp.items())
    ]

    return residues, seen_asym


def get_all_entity_chain_pairs(block) -> list[tuple[str, str]]:
    """Return all polypeptide (entity_id, chain_id) pairs present in the block.

    Entity IDs are collected from ``_entity_poly_seq`` (or ``_entity_poly``
    when that category is absent) and filtered to polypeptides only via
    ``_entity_poly.type``.  For each entity the corresponding author chain IDs
    are read from ``_atom_site`` ATOM records.

    Args:
        block: A :class:`gemmi.cif.Block` object.

    Returns:
        Sorted list of ``(entity_id, chain_id)`` string tuples.
    """
    # Collect unique entity_ids from _entity_poly_seq or _entity_poly
    eps = block.get_mmcif_category("_entity_poly_seq")
    if eps:
        entity_ids = sorted(set(eps.get("entity_id", [])))
    else:
        ep = block.get_mmcif_category("_entity_poly")
        entity_ids = sorted(set(ep.get("entity_id", []))) if ep else []

    # Filter to polypeptides
    ep = block.get_mmcif_category("_entity_poly")
    polypeptide_entities: set[str] = set()
    if ep:
        for idx, eid in enumerate(ep.get("entity_id", [])):
            poly_type = ep.get("type", [])[idx] if ep.get("type") else ""
            if "peptide" in str(poly_type).lower():
                polypeptide_entities.add(eid)
    else:
        polypeptide_entities = set(entity_ids)

    entity_ids = [e for e in entity_ids if e in polypeptide_entities]

    # For each entity, collect unique auth chain IDs from _atom_site
    atom_site = block.get_mmcif_category("_atom_site")
    if not atom_site:
        return []

    groups = atom_site.get("group_PDB", [])
    at_entity_ids = atom_site.get("label_entity_id", [])
    chain_ids_col = atom_site.get("auth_asym_id", [])
    model_nums = atom_site.get("pdbx_PDB_model_num", [])

    entity_to_chains: dict[str, list[str]] = {e: [] for e in entity_ids}
    seen_pairs: set[tuple[str, str]] = set()

    for i, group in enumerate(groups):
        if group != "ATOM":
            continue
        if model_nums and model_nums[i] not in ("1", 1):
            continue
        eid = at_entity_ids[i]
        if eid not in entity_to_chains:
            continue
        pair = (eid, chain_ids_col[i])
        if pair not in seen_pairs:
            seen_pairs.add(pair)
            entity_to_chains[eid].append(chain_ids_col[i])

    pairs: list[tuple[str, str]] = []
    for eid in entity_ids:
        for chain in entity_to_chains[eid]:
            pairs.append((eid, chain))

    return sorted(pairs)


def build_pdbx_poly_seq_scheme(
    block,
    entity_id: str,
    chain_id: str,
    cc: ChemCompMapping,
) -> dict:
    """Reconstruct ``_pdbx_poly_seq_scheme`` for a single (entity, chain) pair.

    Uses a local pairwise alignment (``lalign36``) between the canonical
    deposited sequence and the coordinate sequence to establish the
    correspondence between entity sequence positions and author residue
    numbering.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        entity_id: Entity identifier string (e.g. ``"1"``).
        chain_id: Author chain identifier (e.g. ``"A"``).
        cc: :class:`~pdbe_sifts.mmcif.chem_comp.ChemCompMapping` instance for
            three-letter → one-letter conversion.

    Returns:
        A dict mapping each of the 12 ``_pdbx_poly_seq_scheme`` column names
        to a parallel list of values.  Returns an empty dict when data is
        insufficient.
    """
    canonical = get_canonical_residues(block, entity_id)
    coords, asym_id = get_coordinate_residues(block, entity_id, chain_id)

    if not canonical:
        logger.warning(
            f"[{entity_id}] No sequence data — cannot build poly_seq_scheme"
        )
        return {}

    canonical_1L = "".join(cc.get(mon) for _, mon in canonical)
    coord_1L = "".join(cc.get(r["comp_id"]) for r in coords)

    # Alignment
    best = None
    if coord_1L:
        try:
            alns = list(do_alignment_lalign36(canonical_1L, coord_1L))
            best = alns[0] if alns else None
        except Exception as exc:
            logger.error(f"lalign36 failed: {exc}")

    # Build mapping: canonical index (0-based) → coord dict or None
    n_can = len(canonical)
    mapping: list[dict | None] = [None] * n_can

    if best is not None:
        al_can = str(best[0].seq)
        al_coord = str(best[1].seq)
        can_offset = best[0]._al_start - 1  # 1-based → 0-based
        coord_offset = best[1]._al_start - 1

        can_idx = can_offset
        coord_idx = coord_offset

        for can_ch, crd_ch in zip(al_can, al_coord, strict=False):
            if can_ch != "-" and crd_ch != "-":
                if can_idx < n_can and coord_idx < len(coords):
                    mapping[can_idx] = coords[coord_idx]
                can_idx += 1
                coord_idx += 1
            elif can_ch != "-":
                can_idx += 1
            else:
                coord_idx += 1

    # Build parallel column lists
    cols: dict[str, list] = {k: [] for k in SCHEME_COLUMNS}

    for i, (seq_id, mon_id) in enumerate(canonical):
        obs = mapping[i]
        cols["asym_id"].append(obs["asym_id"] if obs else (asym_id or _MISSING))
        cols["entity_id"].append(entity_id)
        cols["seq_id"].append(str(seq_id))
        cols["mon_id"].append(mon_id)
        cols["ndb_seq_num"].append(str(seq_id))
        cols["pdb_seq_num"].append(str(obs["auth_seq_id"]) if obs else _MISSING)
        cols["auth_seq_num"].append(
            str(obs["auth_seq_id"]) if obs else _MISSING
        )
        cols["pdb_mon_id"].append(obs["comp_id"] if obs else mon_id)
        cols["auth_mon_id"].append(obs["comp_id"] if obs else mon_id)
        cols["pdb_strand_id"].append(chain_id)
        cols["pdb_ins_code"].append(obs["ins_code"] if obs else ".")
        cols["hetero"].append(obs["hetero"] if obs else "n")

    logger.debug(
        f"[{entity_id}/{chain_id}] poly_seq_scheme: "
        f"{sum(v != _MISSING for v in cols['pdb_seq_num'])}/{n_can} observed"
    )
    return cols


def build_pdbx_poly_seq_scheme_all(block, cc: ChemCompMapping) -> dict:
    """Reconstruct ``_pdbx_poly_seq_scheme`` for all polypeptide chains.

    Enumerates every (entity_id, chain_id) pair via
    :func:`get_all_entity_chain_pairs`, calls
    :func:`build_pdbx_poly_seq_scheme` for each, and merges the results into
    a single dict of parallel lists.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        cc: :class:`~pdbe_sifts.mmcif.chem_comp.ChemCompMapping` instance.

    Returns:
        Merged dict with the same 12-column structure as
        :func:`build_pdbx_poly_seq_scheme`.  Returns an empty dict when no
        polypeptide chains are found.
    """
    pairs = get_all_entity_chain_pairs(block)
    if not pairs:
        logger.warning("No polypeptide entity/chain pairs found in block")
        return {}

    logger.info(
        f"Reconstructing _pdbx_poly_seq_scheme for {len(pairs)} pair(s): {pairs}"
    )
    merged: dict[str, list] = {k: [] for k in SCHEME_COLUMNS}

    for entity_id, chain_id in pairs:
        partial = build_pdbx_poly_seq_scheme(block, entity_id, chain_id, cc)
        for col in SCHEME_COLUMNS:
            merged[col].extend(partial.get(col, []))

    return merged


def write_pdbx_poly_seq_scheme(block, scheme: dict) -> None:
    """Write *scheme* into *block* as ``_pdbx_poly_seq_scheme``.

    Args:
        block: A :class:`gemmi.cif.Block` object (modified in place).
        scheme: Dict returned by :func:`build_pdbx_poly_seq_scheme` or
            :func:`build_pdbx_poly_seq_scheme_all`.
    """
    block.set_mmcif_category("_pdbx_poly_seq_scheme", scheme)
