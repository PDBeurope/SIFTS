#!/usr/bin/env python3


import ast
import csv
import gzip
import os
from pathlib import Path
from typing import NamedTuple

from funcy.debug import log_durations

import pdbe_sifts.segments_generation.alignment as align
from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.chain import Chain
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.unp.unp import UNP


class XRefSegment(NamedTuple):
    entry_id: str  # VARCHAR(4 BYTE) NOT NULL ENABLE,
    entity_id: int  # INTEGER NOT NULL ENABLE,
    id: int  # noqa: A003  # INTEGER NOT NULL ENABLE,
    auth_asym_id: str  # VARCHAR(10 BYTE),
    struct_asym_id: str  # VARCHAR(4 BYTE) NOT NULL ENABLE,
    accession: str | None  # VARCHAR(20 BYTE),
    name: str | None  # VARCHAR(50 BYTE),
    seq_version: int | None  # INTEGER,
    unp_start: int | None  # INTEGER,
    pdb_start: int  # INTEGER,
    unp_end: int | None  # INTEGER,
    pdb_end: int  # INTEGER,
    auth_start: int  # INTEGER,
    auth_start_icode: str  # VARCHAR(1 BYTE),
    auth_end: int  # INTEGER,
    auth_end_icode: str  # VARCHAR(1 BYTE),
    conflicts: int | None  # INTEGER,
    modifications: int | None  # INTEGER,
    unp_alignment: str | None  # CLOB,
    pdb_alignment: str  # CLOB,
    identity: float | None  # BINARY_FLOAT,
    score: float | None  # BINARY_FLOAT,
    best_mapping: bool  # BOOLEAN NOT NULL CHECK (BEST_MAPPING in (1,0)),
    canonical_acc: bool  # BOOLEAN NOT NULL CHECK (CANONICAL_ACC in (1,0)),
    reference_acc: str | None  # VARCHAR(20 BYTE),
    chimera: bool  # BOOLEAN


class XRefResidue(NamedTuple):
    entry_id: str  # VARCHAR(4 BYTE) NOT NULL ENABLE,
    entity_id: int  # INTEGER NOT NULL ENABLE,
    id: int  # noqa: A003  # INTEGER NOT NULL ENABLE,
    auth_asym_id: str  # VARCHAR(10 BYTE),
    struct_asym_id: str  # VARCHAR(4 BYTE) NOT NULL ENABLE,
    unp_segment_id: int  # INTEGER NOT NULL ENABLE,
    auth_seq_id: int | None  # INTEGER,
    auth_seq_id_ins_code: str  # VARCHAR(1 BYTE),
    pdb_seq_id: int  # INTEGER,
    unp_seq_id: int  # INTEGER,
    observed: str  # VARCHAR(1 BYTE),
    dbentry_id: int | None  # INTEGER,
    accession: str  # VARCHAR(20 BYTE),
    name: str  # VARCHAR(50 BYTE),
    type: str  # noqa: A003  # VARCHAR(50 BYTE),
    unp_one_letter_code: str  # VARCHAR(1 BYTE),
    pdb_one_letter_code: str  # VARCHAR(4 BYTE),
    chem_comp_id: str  # VARCHAR(3 BYTE),
    mh_id: int  # INTEGER,
    tax_id: int  # INTEGER,
    canonical_acc: bool  # BOOLEAN NOT NULL CHECK (CANONICAL_ACC in (1,0)),
    reference_acc: str  # VARCHAR(20 BYTE),
    best_mapping: bool  # BOOLEAN,
    residue_id: str  # VARCHAR(15 BYTE) NOT NULL ENABLE


@log_durations(logger.debug)
def process_seg_csv_line(row, sifts_seg_col=26):
    """Reads the sifts_xref_segment.csv line and filters the following @params
    pdb, entity, struct_asym_id, accession, unp_start, pdb_start, unp_end, pdb_end
    @input :row, where row is csv line split by "," whose total length must be equal to
    total columns in sifts_xref_segment.csv
    """

    if len(row) == sifts_seg_col:
        (
            pdb,
            entity,
            struct_asym_id,
            accession,
            unp_start,
            pdb_start,
            unp_end,
            pdb_end,
        ) = (row[0], row[1], row[4], row[5], row[8], row[9], row[10], row[11])
        return (
            pdb,
            entity,
            struct_asym_id,
            accession,
            unp_start,
            pdb_start,
            unp_end,
            pdb_end,
        )
    else:
        logger.error(
            "Unexpected number of columns in seg_csv. "
            f"{len(row)} columns instead of {sifts_seg_col}"
        )
        return False


@log_durations(logger.debug)
def read_prev_csv(file_name):
    """
    Reads previous week seg csv and stores following info to compare it to the new data
    returns dict_of_dict storing following @params
    pdb, entity, struct_asym_id, accession, unp_start, pdb_start, unp_end, pdb_end
    """
    data = []
    if os.path.exists(file_name):
        with gzip.open(file_name, "rt") as csvfile:
            for row in csv.reader(csvfile):
                data.append(XRefSegment(*row))
    return data


@log_durations(logger.debug)
def track_mapping_change(data, seg_file_name, res_file_name):
    """
    Compare the new data with prev week's csv and returns are list as output
    @output= track_data, which is a list_of_list with following format
        [pdb, entity, struct_asym_id, accession, unp_start, pdb_start,
        unp_end, pdb_end, mapping_changed_flag, what_changed]
    """

    track_data = []  # list of list
    if os.path.exists(seg_file_name) and os.path.exists(res_file_name):
        prev_info = read_prev_csv(seg_file_name)
    else:
        prev_info = {}
    for row in data:
        delta_flag = False
        delta = ""
        (
            pdb,
            entity,
            struct_asym_id,
            accession,
            unp_start,
            pdb_start,
            unp_end,
            pdb_end,
        ) = process_seg_csv_line(row)
        pdb_start = str(pdb_start)
        if accession:
            if pdb in prev_info:
                if entity in prev_info[pdb]:
                    if struct_asym_id in prev_info[pdb][entity]:
                        if accession in prev_info[pdb][entity][struct_asym_id]:
                            if (
                                pdb_start
                                in prev_info[pdb][entity][struct_asym_id][
                                    accession
                                ]
                            ):
                                prev_seg_list = prev_info[pdb][entity][
                                    struct_asym_id
                                ][accession][pdb_start]
                                new_seg = ",".join(
                                    [str(pdb_end), str(unp_start), str(unp_end)]
                                )
                                if new_seg in prev_seg_list:
                                    # same segments
                                    pass
                                else:
                                    delta_flag = True
                                    delta = "New segment"
                            else:
                                delta_flag = True
                                delta = "New pdb_start"
                        else:
                            logger.info(
                                f"UNP Accession changed, new accession {accession}"
                            )
                            delta_flag = True
                            delta = "New unp accession"

                    else:
                        delta_flag = True
                        delta = "New mapped struct_asym_id"
                else:
                    delta_flag = True
                    delta = "New mapped pdb entity"

            else:
                delta_flag = True
                delta = "New mapped pdb entry"

            if delta_flag:
                track_data.append(
                    (
                        pdb,
                        entity,
                        struct_asym_id,
                        accession,
                        unp_start,
                        pdb_start,
                        unp_end,
                        pdb_end,
                        delta_flag,
                        delta,
                    )
                )

    return track_data


@log_durations(logger.debug)
def write_seg_csv(out_dir, pdbid, data, nf90_mode):
    write_flag = False
    file_name = Path(out_dir, f"{pdbid}_seg.csv.gz")
    nf90_file_name = Path(out_dir, f"{pdbid}_nf90_seg.csv.gz")
    del_file_name = Path(out_dir, f"{pdbid}_delta.csv.gz")
    if nf90_mode:
        logger.info(f"Writing csv for {nf90_file_name}")
        write_flag = True

        with gzip.open(nf90_file_name, "wt") as f:
            sw = csv.writer(
                f,
                delimiter=",",
                quotechar="|",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            for item in data:
                item = fix_values(item)
                sw.writerow(item)
    else:
        # track if uniprot mapping has changed or not
        res_file_name = Path(out_dir, f"{pdbid}_res.csv.gz")
        unp_delta = []
        unp_delta = track_mapping_change(data, file_name, res_file_name)

        # write the csv only when mapping is changed
        if not unp_delta:
            if data:
                logger.info(
                    f"Same mapping! Skipped writing csv for {file_name}"
                )

            # remove the del_file_name from previous week
            if os.path.exists(del_file_name):
                os.remove(del_file_name)
        else:
            logger.info(f"Writing csv for {file_name}")
            write_flag = True
            with gzip.open(file_name, "wt") as f:
                sw = csv.writer(
                    f,
                    delimiter=",",
                    quotechar="|",
                    quoting=csv.QUOTE_MINIMAL,
                    lineterminator="\n",
                )
                for item in data:
                    item = fix_values(item)
                    sw.writerow(item)

            logger.info(f"Writing Uniprot mapping changes {del_file_name}")
            with gzip.open(del_file_name, "wt") as df:
                csv_w = csv.writer(
                    df,
                    delimiter=",",
                    quotechar="|",
                    quoting=csv.QUOTE_MINIMAL,
                    lineterminator="\n",
                )
                for item in unp_delta:
                    csv_w.writerow(item)

    return write_flag


def fix_values(row):
    # Replace commas with ; to not have to quote strings
    row = [
        col.replace(",", ";") if isinstance(col, str) else col for col in row
    ]
    # replace True with 0, False with 1
    row = [int(col) if isinstance(col, bool) else col for col in row]
    return row


@log_durations(logger.debug)
def write_res_csv(out_dir, pdbid, data):
    file_name = Path(out_dir, f"{pdbid}_res.csv.gz")
    logger.info(f"Writing csv for {file_name}")
    with gzip.open(file_name, "wt") as f:
        sw = csv.writer(
            f,
            delimiter=",",
            quotechar="|",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )

        for item in data:
            for item2 in item:
                sw.writerow(fix_values(item2))


@log_durations(logger.debug)
def insert_residues(
    chain_obj: Chain, unp_object=None, iso_acc=None, canonical=False, conn=None
):
    # 23 columns for residues

    if not chain_obj.is_chimera:
        iso = iso_acc
        unp_obj = unp_object

        residue_map = chain_obj.residue_maps.get(iso, {})

        if chain_obj.tax_id is not None and not chain_obj.is_chimera:
            tax_id = chain_obj.tax_id
        elif unp_obj and unp_obj.taxonomy:
            tax_id = (
                int(unp_obj.taxonomy[0]) if unp_obj.taxonomy != [] else None
            )
        else:
            tax_id = None

    res: list[tuple] = []
    index = 1

    for r_list in chain_obj.residues:
        # For MH
        if not isinstance(r_list, list):
            r_list = [r_list]

        for r in r_list:
            if chain_obj.is_chimera:
                residue_map = {}
                mapped = False
                unp_obj = None
                iso = None
                tax_id = None

                for i, u in zip(iso_acc, unp_object, strict=False):
                    if r.n in chain_obj.residue_maps[i]:
                        residue_map = chain_obj.residue_maps[i]
                        mapped = True
                        unp_obj = u
                        iso = i
                        tax_id = (
                            int(unp_obj.taxonomy[0])
                            if unp_obj.taxonomy != []
                            else None
                        )
                        break
            else:
                mapped = r.n in residue_map

            # For Chromophores which map to UniProt we "unroll" the r.oneL
            if mapped and isinstance(residue_map[r.n], list):
                res = process_chromophores(
                    chain_obj,
                    canonical,
                    iso,
                    unp_obj,
                    residue_map,
                    tax_id,
                    res,
                    index,
                    r,
                    mapped,
                )
            else:
                best_flag = False
                if mapped:
                    unp_residue = unp_obj.seq_isoforms[iso][
                        residue_map[r.n] - 1
                    ]

                    # Remove old conflicts which don't exist anymore
                    # (due to a PDB or UniProt sequence update)
                    if r.rtype == "Conflict" and r.oneL == unp_residue:
                        r.rtype = None
                    elif r.rtype is None and r.oneL != unp_residue:
                        # Make a conflict if they are different
                        r.rtype = "Conflict"

                    # best mapping
                    best_flag = (
                        chain_obj.best[chain_obj.canonicals[0]][0] == iso
                        if not chain_obj.is_chimera
                        else True
                    )

                res_id = [
                    chain_obj.pdbid,
                    chain_obj.entity_id,
                    chain_obj.struct_asym_id,
                    r.n,
                ]
                res_id = [str(item) for item in res_id]
                res_id = "_".join(res_id)
                res.append(
                    XRefResidue(
                        chain_obj.pdbid,
                        chain_obj.entity_id,
                        index,
                        chain_obj.auth_asym_id,
                        chain_obj.struct_asym_id,
                        1,  # TODO: UNP_SEGMENT_ID
                        int(r.auth_n) if r.auth_n else None,
                        r.auth_ins if r.auth_ins else " ",
                        r.n,
                        residue_map[r.n] if mapped else None,  # UNP_SEQ_ID
                        "Y" if r.observed else "N",
                        None,  # Future DB ENTRY
                        iso if mapped else None,
                        (
                            unp_obj.longName
                            if unp_obj is not None and mapped
                            else None
                        ),  # NAME
                        r.rtype,  # TYPE
                        unp_residue if mapped else None,  # UNP_ONE_LETTER_CODE
                        r.oneL,
                        r.threeL,
                        r.mh,
                        tax_id if mapped or unp_obj is None else None,  # TAX_ID
                        (
                            iso in chain_obj.canonicals
                            if iso is not None
                            else canonical
                        ),  # CANONICAL
                        iso,  # REFERENCE_ACC
                        best_flag,  # best mapping
                        res_id,
                    )
                )

            index += 1

    return res


@log_durations(logger.debug)
def process_chromophores(
    chain_obj,
    canonical,
    iso,
    unp_obj: UNP,
    residue_map,
    tax_id,
    res,
    index,
    r,
    mapped,
):
    for idx, m in enumerate(residue_map[r.n]):
        unp_residue = unp_obj.seq_isoforms[iso][m - 1]

        # Remove old conflicts which don't exist anymore
        # (due to a PDB or UniProt sequence update)
        if (
            r.rtype == "Conflict"
            and idx <= (len(r.oneL) - 1)
            and r.oneL[idx] == unp_residue
        ):
            r.rtype = None
            # Make a conflict if they are different
        elif (
            r.rtype is None
            and idx <= (len(r.oneL) - 1)
            and r.oneL[idx] != unp_residue
        ):
            r.rtype = "Conflict"
        best_flag = (
            chain_obj.best[chain_obj.canonicals[0]][0] == iso
            if not chain_obj.is_chimera
            else True
        )  # best mapping

        res_id = [
            chain_obj.pdbid,
            chain_obj.entity_id,
            chain_obj.struct_asym_id,
            r.n,
        ]
        res_id = [str(item) for item in res_id]
        res_id = "_".join(res_id)
        res.append(
            XRefResidue(
                chain_obj.pdbid,
                chain_obj.entity_id,
                index,
                chain_obj.auth_asym_id,
                chain_obj.struct_asym_id,
                1,  # TODO: UNP_SEGMENT_ID
                int(r.auth_n) if r.auth_n else None,
                r.auth_ins if r.auth_ins else " ",
                r.n,
                m if mapped else None,  # UNP_SEQ_ID
                "Y" if r.observed else "N",
                None,  # DBENTRY_ID
                iso if mapped else None,
                unp_obj.longName
                if unp_obj is not None and mapped
                else None,  # NAME
                r.rtype,  # TYPE
                unp_residue if mapped else None,  # UNP_ONE_LETTER_CODE
                r.oneL[idx] if idx <= (len(r.oneL) - 1) else r.oneL[0],
                r.threeL,
                r.mh,
                tax_id if mapped or unp_obj is None else None,  # TAX_ID
                (
                    iso in chain_obj.canonicals
                    if iso is not None
                    else canonical
                ),  # CANONICAL
                iso,  # REFERENCE_ACC
                best_flag,
                res_id,
            )
        )
    return res


@log_durations(logger.debug)
def insert_mappings(out_dir, entry_obj: Entry, nf90_mode, conn=None):
    segs = []
    xref_residue = []
    pdbid = entry_obj.pdbid
    for _, chain_obj in list(entry_obj.chains.items()):
        if not nf90_mode and chain_obj.skip:
            logger.warning(
                f"Told to skip chain {chain_obj.auth_asym_id}. Skipping"
            )
            continue
        logger.info(chain_obj.mappings)
        # If there are no segments OR the canonical couldn't be mapped
        # (e.g. 1loi).
        # We have to ensure a CANONICAL_ACC = 1 for each chain so,
        # if the canonical didn't align, we create an "empty" mapping
        if not chain_obj.segments or (
            not any(x in chain_obj.segments for x in chain_obj.canonicals)
        ):
            auth_s, ins_s = chain_obj.get_residue_auth(1)
            auth_e, ins_e = chain_obj.get_residue_auth(len(chain_obj.sequence))

            segs.append(
                XRefSegment(
                    chain_obj.pdbid,
                    chain_obj.entity_id,
                    1,
                    chain_obj.auth_asym_id,
                    chain_obj.struct_asym_id,
                    None,
                    None,
                    None,
                    None,
                    1,
                    None,
                    len(chain_obj.sequence),
                    auth_s,
                    ins_s,
                    auth_e,
                    ins_e,
                    None,
                    None,
                    None,
                    chain_obj.sequence,
                    None,
                    None,
                    len(chain_obj.segments) == 0,
                    True,
                    None,
                    chain_obj.is_chimera,
                )
            )

            if not nf90_mode:
                xref_residue.append(
                    insert_residues(chain_obj, canonical=True, conn=conn)
                )

        # A shared one for all the accessions
        if chain_obj.is_chimera:
            blanks = [(1, len(chain_obj.sequence))]

        # Each mapping represents a different UniProt isoform
        for iso, maps in list(chain_obj.segments.items()):
            col_id = 1
            # A set of blank segments for each accession
            if not chain_obj.is_chimera:
                blanks = [(1, len(chain_obj.sequence))]

            unp_obj = entry_obj.accessions[iso]

            # get the alignment object for that isoform
            al = chain_obj.mappings[iso][0][2]

            blanks, col_id = process_non_chimera(
                nf90_mode,
                segs,
                xref_residue,
                chain_obj,
                iso,
                maps,
                col_id,
                unp_obj,
                al,
                blanks,
                conn,
            )

        if chain_obj.is_chimera:
            handle_chimera(segs, chain_obj, blanks, iso, col_id)

            iso = list(chain_obj.segments.keys())
            unp_obj = [entry_obj.accessions[i] for i in iso]

            if not nf90_mode:
                xref_residue.append(
                    insert_residues(
                        chain_obj,
                        unp_object=unp_obj,
                        iso_acc=iso,
                        canonical=True,
                        conn=conn,
                    )
                )

    # only do this if there are segments to insert
    write_flag = False

    write_flag = write_seg_csv(out_dir, pdbid, segs, nf90_mode)
    if write_flag:
        if xref_residue:
            write_res_csv(out_dir, pdbid, xref_residue)
    else:
        logger.info(f"Same mapping! Skipped writing res csv for {pdbid}")

    return segs, xref_residue


@log_durations(logger.debug)
def process_non_chimera(
    nf90_mode,
    segs,
    xref_residue,
    chain_obj,
    iso,
    maps,
    col_id,
    unp_obj,
    al,
    blanks,
    conn,
):
    for pdb_range, unp_range in maps:
        logger.info(
            f"Insert mappings: [{chain_obj.entity_id}] "
            f"[{chain_obj.auth_asym_id}] [{iso}]: {pdb_range} -> {unp_range}"
        )

        unp_seq = align.get_align_chunk(
            al[0].seq,
            unp_range[0] - al[0]._al_start + 1,
            unp_range[1] - al[0]._al_start + 1,
        )
        pdb_seq = align.get_align_chunk(
            al[1].seq,
            pdb_range[0] - al[1]._al_start + 1,
            pdb_range[1] - al[1]._al_start + 1,
        )

        auth_s, ins_s = chain_obj.get_residue_auth(pdb_range[0])
        auth_e, ins_e = chain_obj.get_residue_auth(pdb_range[1])

        # PDBE-3790, incorrect seq identity esp. in case of same acc repeats
        my_score = chain_obj.scores[iso][0]
        my_identity = chain_obj.scores[iso][1]

        for my_seg in chain_obj.seg_scores[iso]:
            my_range = ast.literal_eval(my_seg)[0]
            # print pdb_range
            if my_range[0] <= pdb_range[0] and pdb_range[1] <= my_range[1]:
                my_score = chain_obj.seg_scores[iso][my_seg]
                my_identity = my_score

        segs.append(
            XRefSegment(
                chain_obj.pdbid,
                chain_obj.entity_id,
                col_id,
                chain_obj.auth_asym_id,
                chain_obj.struct_asym_id,
                iso,
                unp_obj.longName,
                unp_obj.date_seq_update[1],
                unp_range[0],
                pdb_range[0],
                unp_range[1],
                pdb_range[1],
                auth_s,
                ins_s,
                auth_e,
                ins_e,
                align.get_conflicts(unp_seq, pdb_seq),
                None,
                unp_seq,
                pdb_seq,
                my_score,
                my_identity,
                (
                    chain_obj.best[chain_obj.canonicals[0]][0] == iso
                    if not chain_obj.is_chimera
                    else True
                ),
                iso in chain_obj.canonicals,
                iso,
                chain_obj.is_chimera,
            )
        )

        col_id += 1
        blanks, _ = align.remove_range_alignment(
            blanks, blanks, pdb_range, "", False
        )

    if not chain_obj.is_chimera:
        for b in blanks:
            pdb_seq = chain_obj.sequence[b[0] - 1 : b[1]]

            auth_s, ins_s = chain_obj.get_residue_auth(b[0])
            auth_e, ins_e = chain_obj.get_residue_auth(b[1])

            segs.append(
                XRefSegment(
                    chain_obj.pdbid,
                    chain_obj.entity_id,
                    col_id,
                    chain_obj.auth_asym_id,
                    chain_obj.struct_asym_id,
                    None,
                    None,
                    None,
                    None,
                    b[0],
                    None,
                    b[1],
                    auth_s,
                    ins_s,
                    auth_e,
                    ins_e,
                    None,
                    None,
                    None,
                    pdb_seq,
                    None,
                    None,
                    chain_obj.best[chain_obj.canonicals[0]][0] == iso,
                    iso in chain_obj.canonicals,
                    iso,
                    chain_obj.is_chimera,
                )
            )

            col_id += 1

        if not nf90_mode:
            xref_residue.append(
                insert_residues(
                    chain_obj, unp_object=unp_obj, iso_acc=iso, conn=conn
                )
            )

    return blanks, col_id


@log_durations(logger.debug)
def handle_chimera(segs, chain_obj, blanks, iso, col_id):
    for b in blanks:
        pdb_seq = chain_obj.sequence[b[0] - 1 : b[1]]

        auth_s, ins_s = chain_obj.get_residue_auth(b[0])
        auth_e, ins_e = chain_obj.get_residue_auth(b[1])

        segs.append(
            XRefSegment(
                chain_obj.pdbid,
                chain_obj.entity_id,
                col_id,
                chain_obj.auth_asym_id,
                chain_obj.struct_asym_id,
                None,
                None,
                None,
                None,
                b[0],
                None,
                b[1],
                auth_s,
                ins_s,
                auth_e,
                ins_e,
                None,
                None,
                None,
                pdb_seq,
                None,
                None,
                True,
                True,
                iso,
                chain_obj.is_chimera,
            )
        )

        col_id += 1
