import csv
import gzip
import os

from pdbe_sifts.base.log import logger


def get_unp_segments(pdbid, seg_csv, pdbecursor):
    res_info = {}

    sifts_seg_header = [
        "entry_id",
        "entity_id",
        "id",
        "auth_asym_id",
        "struct_asym_id",
        "accession",
        "name",
        "seq_version",
        "unp_start",
        "pdb_start",
        "unp_end",
        "pdb_end",
        "auth_start",
        "auth_start_icode",
        "auth_end",
        "auth_end_icode",
        "conflicts",
        "modifications",
        "unp_alignment",
        "pdb_alignment",
        "identity",
        "score",
        "best_mapping",
        "canonical_acc",
        "reference_acc",
        "chimera",
    ]

    data = {}
    reader = []

    if seg_csv:
        if not os.path.exists(seg_csv):
            logger.warning(f"No sifts seg csv found for {pdbid}")
        else:
            f = gzip.open(seg_csv, "rt")
            csvfile = f.readlines()
            f.close()
            reader = csv.DictReader(csvfile, delimiter=",", fieldnames=sifts_seg_header)
    else:
        rows = pdbecursor.execute(
            "select * from sifts_xref_segment where entry_id=?", [pdbid]
        ).fetchdf()
        reader = [row[1].to_dict() for row in rows.iterrows()]

    for row in reader:
        if row["accession"]:
            data.setdefault(
                "_".join(
                    [
                        row["entry_id"],
                        str(row["entity_id"]),
                        row["struct_asym_id"],
                        row["accession"],
                    ]
                ),
                {},
            ).setdefault(f"{row['unp_start']}_{row['unp_end']}", []).append(
                "_".join(
                    [
                        str(row["pdb_start"]),
                        str(row["pdb_end"]),
                        '1' if row["best_mapping"] else '0',
                        str(row["identity"]),
                    ]
                )
            )

    cif_entity_id = []
    cif_asym_id = []
    cif_unp_accession = []
    cif_segment_id = []
    cif_instance_id = []
    cif_unp_start = []
    cif_unp_end = []
    cif_seq_id_start = []
    cif_seq_id_end = []
    cif_best_mapping = []
    cif_identity = []
    # finding seg_id,insnce_id for a given pdb_asymid_unpacc
    for row in sorted(data):
        seg_boo, ins_boo = 1, 1
        for seg in data[row]:
            for copy in data[row][seg]:
                pdbid, entity_id, asym_id, unp_acc = row.split("_")
                unp_start, unp_end = seg.split("_")
                pdb_start, pdb_end, best_mapping, identity = copy.split("_")
                identity = round(float(identity), 3)
                best_mapping = int(best_mapping)
                best_mapping = "y" if best_mapping else "n"

                cif_entity_id.append(entity_id)
                cif_asym_id.append(asym_id)
                cif_unp_accession.append(unp_acc)
                cif_segment_id.append(seg_boo)
                cif_instance_id.append(ins_boo)
                cif_unp_start.append(unp_start)
                cif_unp_end.append(unp_end)
                cif_seq_id_start.append(pdb_start)
                cif_seq_id_end.append(pdb_end)
                cif_best_mapping.append(best_mapping)
                cif_identity.append(identity)
                res_info.setdefault(int(entity_id), {}).setdefault(
                    asym_id, {}
                ).setdefault(int(pdb_start), []).append(
                    (int(pdb_end), unp_acc, seg_boo, ins_boo)
                )

                if len(data[row][seg]) != 1:
                    ins_boo = ins_boo + 1
                else:
                    seg_boo = seg_boo + 1

    mmcif_cat = (
        cif_entity_id,
        cif_asym_id,
        cif_unp_accession,
        cif_segment_id,
        cif_instance_id,
        cif_unp_start,
        cif_unp_end,
        cif_seq_id_start,
        cif_seq_id_end,
        cif_best_mapping,
        cif_identity,
    )

    return mmcif_cat, res_info


def get_unpres_mapping(pdbid, res_csv, pdbecursor):
    sifts_res_header = [
        "entry_id",
        "entity_id",
        "id",
        "auth_asym_id",
        "struct_asym_id",
        "unp_segment_id",
        "auth_seq_id",
        "auth_seq_id_ins",
        "pdb_seq_id",
        "unp_seq_id",
        "observed",
        "dbentry_id",
        "accession",
        "name",
        "type",
        "unp_one_letter_code",
        "pdb_one_letter_code",
        "chem_comp_id",
        "mh_id",
        "tax_id",
        "canonical_acc",
        "reference_acc",
        "best_mapping",
    ]

    data = {}
    mon_id = {}
    reader = []

    if res_csv:
        if not os.path.exists(res_csv):
            logger.warning(f"No sifts res csv found for {pdbid}")
        else:
            f = gzip.open(res_csv, "rt")
            csvfile = f.readlines()
            f.close()
            reader = csv.DictReader(csvfile, delimiter=",", fieldnames=sifts_res_header)
    else:
        rows = pdbecursor.execute(
            "select * from sifts_xref_residue where entry_id=?", [pdbid]
        ).fetchdf()
        reader = [row[1].to_dict() for row in rows.iterrows()]

    for row in reader:
        if row["accession"] and row["best_mapping"]:
            data.setdefault(int(row["entity_id"]), {}).setdefault(
                row["struct_asym_id"], {}
            ).setdefault(int(row["pdb_seq_id"]), []).append(
                (
                    row["chem_comp_id"],
                    row["pdb_one_letter_code"],
                    row["unp_one_letter_code"],
                    row["accession"],
                    int(row["unp_seq_id"]),
                    row["type"],
                    row["mh_id"],
                    row["observed"].lower(),
                )
            )
        mon_id.setdefault(int(row["entity_id"]), {}).setdefault(
            row["struct_asym_id"], {}
        ).setdefault(int(row["pdb_seq_id"]), {})[row["chem_comp_id"]] = row[
            "pdb_one_letter_code"
        ]

    return data, mon_id
