import numpy as np

from pdbe_sifts.base.log import logger


def get_instance(data):
    """
    get instance_id for xref_db
    This is specific to cath/scop2 db as they already have segment_id ie ordinal computed
    """
    for_mmcif = {}  # k=entity_asym_xrefdb_acc, v=list for remaining lines
    for tag in data:
        inst_id = 1
        for domain_name in sorted(data[tag]):
            my_start = sorted(data[tag][domain_name])
            # sort in terms of sequence start id
            for i in range(0, len(my_start)):
                seq_id_start = my_start[i]
                copy_tag = data[tag][domain_name][seq_id_start]
                # print(seq_id_start, copy_tag)
                for copy in copy_tag:
                    seg_id, seq_id_end = copy.split("@")
                    (
                        entity_id,
                        asym_id,
                        xref_db,
                        domain_acc,
                    ) = tag.split("@")
                    entity_id = int(entity_id)
                    # print(entity_id,asym_id,domain_acc,domain_name,seg_id,inst_id)
                    for_mmcif.setdefault(entity_id, {}).setdefault(
                        asym_id, {}
                    ).setdefault(xref_db, []).append(
                        [
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ]
                    )
            # instance increases coz
            # its different domain name but same domain acc
            inst_id = inst_id + 1

    return for_mmcif


def get_cath_info(pdbid, mycursor):
    # getting CATH data
    cath_query = """
        select entity_id,
            struct_asym_id,
            accession,
            domain,
            ordinal,
            "START" as pdb_start,
            "END" as pdb_end
        from
            entity_cath
        where entry_id = :1
    """
    result = mycursor.execute(cath_query, (pdbid,))
    data = {}
    xref_db = "CATH"
    for row in result:
        entity_id = row["entity_id"]
        asym_id = row["struct_asym_id"]
        domain_acc = row["accession"]
        domain_name = row["domain"]
        seg_id = int(row["ordinal"])
        seq_id_start = row["pdb_start"]
        seq_id_end = row["pdb_end"]
        data.setdefault(f"{entity_id}@{asym_id}@{xref_db}@{domain_acc}", {}).setdefault(
            domain_name, {}
        ).setdefault(int(seq_id_start), []).append(f"{seg_id}@{seq_id_end}")

    mmcif_data = get_instance(data)
    return mmcif_data


def get_scop2fa_info(pdbid, mycursor):
    # getting scop2 data
    scop2_query = """
        SELECT
            distinct entity_id, auth_asym_id, struct_asym_id,
            fa_domid , ordinal, pdb_start, pdb_end
        FROM
            entity_scop2_fa
        WHERE
            entry_id= :1
        ORDER BY
            fa_domid,ordinal
        """

    scop2_rows = mycursor.execute(scop2_query, (pdbid,))
    data = {}
    xref_db = "SCOP2"
    for row in scop2_rows:
        entity_id = row["entity_id"]
        asym_id = row["struct_asym_id"]
        domain_acc = row["fa_domid"]
        domain_name = "FA"
        seg_id = int(row["ordinal"])
        seq_id_start = row["pdb_start"]
        seq_id_end = row["pdb_end"]
        data.setdefault(f"{entity_id}@{asym_id}@{xref_db}@{domain_acc}", {}).setdefault(
            domain_name, {}
        ).setdefault(int(seq_id_start), []).append(f"{seg_id}@{seq_id_end}")
    mmcif_data = get_instance(data)
    return mmcif_data


def get_scop2sf_info(pdbid, mycursor):
    # getting scop2 data
    scop2_query = """
    SELECT
        distinct entity_id, auth_asym_id, struct_asym_id, sf_domid , ordinal, pdb_start, pdb_end
    FROM
        entity_scop2_sf
    WHERE
        entry_id= :1
    ORDER BY
        sf_domid,ordinal
    """

    scop2_rows = mycursor.execute(scop2_query, (pdbid,))
    data = {}
    xref_db = "SCOP2"
    for row in scop2_rows:
        logger.debug(row)
        entity_id = row["entity_id"]
        asym_id = row["struct_asym_id"]
        domain_acc = row["sf_domid"]
        domain_name = "SF"
        seg_id = int(row["ordinal"])
        seq_id_start = row["pdb_start"]
        seq_id_end = row["pdb_end"]
        data.setdefault(f"{entity_id}@{asym_id}@{xref_db}@{domain_acc}", {}).setdefault(
            domain_name, {}
        ).setdefault(int(seq_id_start), []).append(f"{seg_id}@{seq_id_end}")
    mmcif_data = get_instance(data)
    return mmcif_data


def get_scop2bsf_info(pdbid, mycursor):
    # getting scop2 data
    scop2_query = """
    SELECT
        distinct entity_id, auth_asym_id, struct_asym_id, sf_domid , ordinal, pdb_start, pdb_end
    FROM
        entity_scop2b_sf
    WHERE
        entry_id= :1
    ORDER BY
        sf_domid,ordinal
    """

    scop2_rows = mycursor.execute(scop2_query, (pdbid,))
    data = {}
    xref_db = "SCOP2B"
    for row in scop2_rows:
        entity_id = row["entity_id"]
        asym_id = row["struct_asym_id"]
        domain_acc = row["sf_domid"]
        domain_name = "SF"
        seg_id = int(row["ordinal"])
        seq_id_start = row["pdb_start"]
        seq_id_end = row["pdb_end"]
        data.setdefault(f"{entity_id}@{asym_id}@{xref_db}@{domain_acc}", {}).setdefault(
            domain_name, {}
        ).setdefault(int(seq_id_start), []).append(f"{seg_id}@{seq_id_end}")
    mmcif_data = get_instance(data)
    return mmcif_data


def get_pfam_info(pdbid, mycursor):
    #  instance id calculation is different from scop2/cath which already have ordinals

    pfam_query = """
        SELECT
            entity_id,
            struct_asym_id,
            accession as accession,
            "START" as pdb_start,
            "END" as pdb_end
        FROM entity_pfam
        WHERE entry_id = :1
    """
    logger.debug(pfam_query)
    pfam_rows = mycursor.execute(pfam_query, (pdbid,))

    xref_db = "Pfam"
    data = {}
    for entity_id, asym_id, domain_acc, seq_id_start, seq_id_end in pfam_rows:
        domain_name = None
        data.setdefault(
            f"{entity_id}@{asym_id}@{xref_db}@{domain_acc}@{domain_name}", {}
        )[seq_id_start] = seq_id_end

    mmcif_data = {}
    # get instance id
    for tag in data:
        entity_id, asym_id, xref_db, domain_acc, domain_name = tag.split("@")
        if domain_name == "None":
            domain_name = None
        entity_id = int(entity_id)
        seg_id, inst_id = 1, 1
        for seq_id_start in sorted(data[tag]):
            seq_id_end = data[tag][seq_id_start]
            mmcif_data.setdefault(entity_id, {}).setdefault(asym_id, {}).setdefault(
                xref_db, []
            ).append(
                [domain_acc, domain_name, seg_id, inst_id, seq_id_start, seq_id_end]
            )
            inst_id = inst_id + 1

    return mmcif_data


def get_xref_info(all_ent, pdbid, pdbecursor):
    cath_data = get_cath_info(pdbid, pdbecursor)
    scop2fa_data = get_scop2fa_info(pdbid, pdbecursor)
    scop2sf_data = get_scop2sf_info(pdbid, pdbecursor)
    scop2bsf_data = get_scop2bsf_info(pdbid, pdbecursor)
    pfam_data = get_pfam_info(pdbid, pdbecursor)

    # sort the data per entity per chain per db
    # db order followed cath, scop2, pfam
    all_data = []
    res_info = {}
    for entity in all_ent:
        for chain in all_ent[entity]:
            # cath db
            if entity in cath_data and chain in cath_data[entity]:
                for db in cath_data[entity][chain]:
                    for tag in cath_data[entity][chain][db]:
                        all_data.append([entity, chain, db] + tag)
                        (
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ) = tag
                        start, end = int(seq_id_start), int(seq_id_end)
                        xref_db = "CATH"
                        res_info.setdefault(entity, {}).setdefault(
                            chain, {}
                        ).setdefault(xref_db, {}).setdefault(start, []).append(
                            (end, domain_acc, domain_name, seg_id, inst_id)
                        )

            # scop2fa
            if entity in scop2fa_data and chain in scop2fa_data[entity]:
                for db in scop2fa_data[entity][chain]:
                    for tag in scop2fa_data[entity][chain][db]:
                        all_data.append([entity, chain, db] + tag)
                        (
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ) = tag
                        start, end = int(seq_id_start), int(seq_id_end)
                        xref_db, domain_name = "SCOP2", "SCOP2-FA"
                        res_info.setdefault(entity, {}).setdefault(
                            chain, {}
                        ).setdefault(xref_db, {}).setdefault(start, []).append(
                            (end, domain_acc, domain_name, seg_id, inst_id)
                        )
            # scop2sf
            if entity in scop2sf_data and chain in scop2sf_data[entity]:
                for db in scop2sf_data[entity][chain]:
                    for tag in scop2sf_data[entity][chain][db]:
                        all_data.append([entity, chain, db] + tag)
                        (
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ) = tag
                        start, end = int(seq_id_start), int(seq_id_end)
                        xref_db, domain_name = "SCOP2", "SCOP2-SF"
                        res_info.setdefault(entity, {}).setdefault(
                            chain, {}
                        ).setdefault(xref_db, {}).setdefault(start, []).append(
                            (end, domain_acc, domain_name, seg_id, inst_id)
                        )
            # scop2bsf
            if entity in scop2bsf_data and chain in scop2bsf_data[entity]:
                for db in scop2bsf_data[entity][chain]:
                    for tag in scop2bsf_data[entity][chain][db]:
                        all_data.append([entity, chain, db] + tag)
                        (
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ) = tag
                        start, end = int(seq_id_start), int(seq_id_end)
                        xref_db, domain_name = "SCOP2", "SCOP2B-SF"
                        res_info.setdefault(entity, {}).setdefault(
                            chain, {}
                        ).setdefault(xref_db, {}).setdefault(start, []).append(
                            (end, domain_acc, domain_name, seg_id, inst_id)
                        )
            # pfam
            if entity in pfam_data and chain in pfam_data[entity]:
                for db in pfam_data[entity][chain]:
                    for tag in pfam_data[entity][chain][db]:
                        all_data.append([entity, chain, db] + tag)
                        (
                            domain_acc,
                            domain_name,
                            seg_id,
                            inst_id,
                            seq_id_start,
                            seq_id_end,
                        ) = tag
                        start, end = int(seq_id_start), int(seq_id_end)
                        xref_db = "Pfam"
                        res_info.setdefault(entity, {}).setdefault(
                            chain, {}
                        ).setdefault(xref_db, {}).setdefault(start, []).append(
                            (end, domain_acc, domain_name, seg_id, inst_id)
                        )

    np_array = np.array(all_data)
    transpose = np_array.T
    transpose_list = transpose.tolist()
    return transpose_list, res_info
