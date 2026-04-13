import argparse
import logging
import re
from collections import defaultdict
from typing import Any

import gemmi

from pdbe_sifts.mmcif.curator.pdbx_connectivity import GapConnectivityChecker
from pdbe_sifts.seq2seq import Seq2Seq

logger = logging.getLogger(__name__)


class CifSequenceUpdater:
    def __init__(self, input_path: str, output_path: str):
        self.input_path = input_path
        self.output_path = output_path
        self.doc: gemmi.cif.Document = gemmi.cif.read(input_path)
        self.block: gemmi.cif.Block = self.doc.sole_block()

    def _pad_alignment(self, alignment_data: dict[str, Any]) -> tuple[str, str]:
        """
        Adjusts alignment sequences to account for missing residues at start/end.
        """
        alignment = alignment_data["alignment"]
        ref_seq_full: str = str(alignment_data["canonical"])

        al_ref: str = str(alignment[0].seq)
        al_coord: str = str(alignment[1].seq)

        start_idx: int = alignment[0]._al_start
        stop_idx: int = alignment[0]._al_stop

        # Prepend missing prefix
        if start_idx > 1:
            num_prefix = start_idx - 1
            al_ref = ref_seq_full[:num_prefix] + al_ref
            al_coord = ("-" * num_prefix) + al_coord

        # Append missing suffix
        if stop_idx < len(ref_seq_full):
            al_ref = al_ref + ref_seq_full[stop_idx:]
            al_coord = al_coord + ("-" * (len(ref_seq_full) - stop_idx))

        return al_ref, al_coord

    def _get_polymer_mapping(self) -> list[dict[str, str]]:
        mappings = []
        entity_poly = self.block.find_mmcif_category("_entity_poly.")
        auth_map = {}

        if entity_poly:
            eid_col = entity_poly.find_column("entity_id")
            strand_col = entity_poly.find_column("pdbx_strand_id")
            for i in range(len(entity_poly)):
                auth_map[eid_col[i]] = [
                    s.strip() for s in strand_col[i].split(",")
                ]

        asym_table = self.block.find_mmcif_category("_struct_asym.")
        if not asym_table:
            return mappings

        ent_id_col = asym_table.find_column("entity_id")
        label_id_col = asym_table.find_column("id")
        used_counts: dict[str, int] = {}

        for i in range(len(asym_table)):
            eid, lid = ent_id_col[i], label_id_col[i]
            if eid in auth_map:
                chains = auth_map[eid]
                idx = used_counts.get(eid, 0)
                if idx < len(chains):
                    mappings.append(
                        {
                            "entity_id": eid,
                            "label_asym_id": lid,
                            "auth_asym_id": chains[idx],
                        }
                    )
                    used_counts[eid] = idx + 1
        return mappings

    def _build_coord_to_auth_seq(
        self, label_asym_id: str, entity_id: str
    ) -> list[tuple[int, str]]:
        table = self.block.find_mmcif_category("_atom_site.")
        if not table:
            return []

        asym_col = table.find_column("label_asym_id")
        entity_col = table.find_column("label_entity_id")
        group_col = table.find_column("group_PDB")
        auth_seq_col = table.find_column("auth_seq_id")
        comp_id_col = table.find_column("label_comp_id")

        coord_to_auth = []
        seen_keys = {}
        coord_idx = 0

        for i in range(len(table)):
            if (
                asym_col[i] != label_asym_id
                or entity_col[i] != entity_id
                or group_col[i] not in ("ATOM", "HETATM")
            ):
                continue
            res_key = (auth_seq_col[i], comp_id_col[i])
            if res_key not in seen_keys:
                coord_idx += 1
                seen_keys[res_key] = coord_idx
                coord_to_auth.append((coord_idx, auth_seq_col[i]))
        return coord_to_auth

    def _parse_raw_seq_mon_ids(self, entity_id: str) -> dict[int, str]:
        """Parse ``_entity_poly.pdbx_seq_one_letter_code`` into ``{seq_id: three_letter}``.

        Standard single-letter residues are expanded via :func:`gemmi.expand_one_letter`.
        Modified residues in parenthetical notation — ``(HYP)``, ``(MSE)``, etc. — are
        extracted directly, preserving the original deposited name.

        Returns an empty dict if the field is absent or the entity is not found.
        """
        cat = self.block.find_mmcif_category("_entity_poly.")
        if not cat:
            return {}
        eid_col = cat.find_column("entity_id")
        raw_col = cat.find_column("pdbx_seq_one_letter_code")
        if not raw_col:
            return {}

        raw_seq = None
        for i in range(len(cat)):
            if eid_col[i] == entity_id:
                raw_seq = raw_col[i]
                break
        if not raw_seq:
            return {}

        # gemmi strips surrounding quotes/semicolons; remove any remaining whitespace
        raw_seq = re.sub(r"\s+", "", raw_seq)

        result: dict[int, str] = {}
        pos = 1
        for m in re.finditer(r"\(([^)]+)\)|([A-Za-z])", raw_seq):
            if m.group(1):  # (HYP) → 'HYP'
                result[pos] = m.group(1).upper()
            else:  # 'P'   → 'PRO'
                result[pos] = gemmi.expand_one_letter(
                    m.group(2).upper(), gemmi.ResidueKind.AA
                )
            pos += 1
        return result

    def _build_alignment_model(
        self, alignment_results: dict[tuple, Any]
    ) -> list[dict[str, Any]]:
        entities = defaultdict(list)
        for (eid, lid, aid), data in alignment_results.items():
            entities[eid].append((lid, aid, data))

        models = []
        for entity_id, chain_list in entities.items():
            _, _, first_data = chain_list[0]
            if "merged_ref" in first_data:
                can_seq_gapped = first_data["merged_ref"]
                coord_seq_gapped = first_data["merged_coord"]
            else:
                can_seq_gapped, coord_seq_gapped = self._pad_alignment(
                    first_data
                )

            coord_to_canonical = {}
            can_idx = coord_idx = 0
            for can_ch, coord_ch in zip(
                can_seq_gapped, coord_seq_gapped, strict=False
            ):
                if can_ch != "-":
                    can_idx += 1
                if coord_ch != "-":
                    coord_idx += 1
                if can_ch != "-" and coord_ch != "-":
                    coord_to_canonical[coord_idx] = can_idx

            canonical_sequence = can_seq_gapped.replace("-", "")
            orig_mon = self._parse_raw_seq_mon_ids(entity_id)
            mon_ids = [
                orig_mon.get(
                    seq_id, gemmi.expand_one_letter(aa, gemmi.ResidueKind.AA)
                )
                for seq_id, aa in enumerate(canonical_sequence, start=1)
            ]

            chains = []
            for lid, aid, _ in chain_list:
                observed_auth = dict(
                    self._build_coord_to_auth_seq(lid, entity_id)
                )
                can_to_coord = {v: k for k, v in coord_to_canonical.items()}

                residues = []
                for seq_id, mon_id in enumerate(mon_ids, start=1):
                    coord_pos = can_to_coord.get(seq_id)
                    auth_seq_id = (
                        observed_auth.get(coord_pos) if coord_pos else None
                    )
                    residues.append(
                        {
                            "seq_id": seq_id,
                            "mon_id": mon_id,
                            "auth_seq_id": auth_seq_id,
                            "observed": auth_seq_id is not None,
                        }
                    )

                chains.append(
                    {
                        "label_asym_id": lid,
                        "auth_asym_id": aid,
                        "residues": residues,
                        "coord_to_canonical": coord_to_canonical,
                    }
                )

            models.append(
                {
                    "entity_id": entity_id,
                    "canonical_sequence": canonical_sequence,
                    "chains": chains,
                }
            )
        return models

    def _write_updates(self, models: list[dict[str, Any]]) -> None:
        # 1. _entity_poly_seq
        self.block.find_mmcif_category("_entity_poly_seq.").erase()
        ep_cat = self.block.init_loop(
            "_entity_poly_seq.", ["entity_id", "num", "mon_id", "hetero"]
        )
        for ent in models:
            for seq_id, mon_id in enumerate(
                [r["mon_id"] for r in ent["chains"][0]["residues"]], start=1
            ):
                ep_cat.add_row([ent["entity_id"], str(seq_id), mon_id, "n"])

        # 2. _pdbx_poly_seq_scheme
        prefix = "_pdbx_poly_seq_scheme."
        tags = [
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
        ]
        self.block.find_mmcif_category(prefix).erase()
        scheme_loop = self.block.init_loop(prefix, tags)
        for ent in models:
            for ch in ent["chains"]:
                for res in ch["residues"]:
                    obs = res["observed"]
                    scheme_loop.add_row(
                        [
                            ch["label_asym_id"],
                            ent["entity_id"],
                            str(res["seq_id"]),
                            res["mon_id"],
                            str(res["seq_id"]),
                            str(res["seq_id"]),
                            res["auth_seq_id"] if obs else ".",
                            res["mon_id"] if obs else ".",
                            res["mon_id"] if obs else ".",
                            ch["auth_asym_id"],
                            ".",
                            "n",
                        ]
                    )

        # 3. _atom_site
        table = self.block.find_mmcif_category("_atom_site.")
        if table:
            asym_col = table.find_column("label_asym_id")
            entity_id_col = table.find_column("label_entity_id")
            seq_id_col = table.find_column("label_seq_id")
            group_col = table.find_column("group_PDB")
            auth_seq_col = table.find_column("auth_seq_id")
            comp_id_col = table.find_column("label_comp_id")

            chain_lookup = {
                c["label_asym_id"]: c["coord_to_canonical"]
                for e in models
                for c in e["chains"]
            }
            asym_to_entity = {
                c["label_asym_id"]: e["entity_id"]
                for e in models
                for c in e["chains"]
            }
            counts, lasts = {}, {}

            for i in range(len(table)):
                lid = asym_col[i]
                if (
                    lid not in chain_lookup
                    or entity_id_col[i] != asym_to_entity.get(lid)
                    or group_col[i] not in ("ATOM", "HETATM")
                ):
                    continue
                key = (auth_seq_col[i], comp_id_col[i])
                if lasts.get(lid) != key:
                    counts[lid] = counts.get(lid, 0) + 1
                    lasts[lid] = key
                c_idx = chain_lookup[lid].get(counts[lid])
                seq_id_col[i] = str(c_idx) if c_idx else "."

    def process(self) -> None:
        """Main execution flow."""
        mappings = self._get_polymer_mapping()
        results = {}
        results_list = []
        for m in mappings:
            data = Seq2Seq(
                self.input_path, m["entity_id"], m["auth_asym_id"]
            ).run()
            data["entity_id"] = m["entity_id"]
            data["chain_id"] = m["auth_asym_id"]
            results[(m["entity_id"], m["label_asym_id"], m["auth_asym_id"])] = (
                data
            )
            results_list.append(data)

        merged = GapConnectivityChecker(
            self.input_path, results_list
        ).get_merged_alignment()

        for (eid, _, aid), data in results.items():
            ma = merged.get((eid, aid))
            if ma:
                data["merged_ref"] = ma["canonical"]
                data["merged_coord"] = ma["coord"]
                if ma["n_fixes"]:
                    logger.info(
                        "Connectivity fix applied: entity=%s chain=%s (%d chunk(s) corrected)",
                        eid,
                        aid,
                        ma["n_fixes"],
                    )

        models = self._build_alignment_model(results)
        self._write_updates(models)
        self.doc.write_file(self.output_path)


def run():
    parser = argparse.ArgumentParser(
        description="Update CIF label_seq_id based on SIFTS alignment."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input CIF file path"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output CIF file path"
    )

    args = parser.parse_args()

    updater = CifSequenceUpdater(args.input, args.output)
    updater.process()


if __name__ == "__main__":
    run()
