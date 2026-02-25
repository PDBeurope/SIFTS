import numpy as np
from gemmi import cif

from pdbe_sifts.base.log import logger


def add_mmcif_cat(cif_cat_obj, my_cat_name, my_cat_val):
    """
    Add mmcif category to mmcif given the [item names],[item values]
    """
    mmcif_cat = cif_cat_obj.setItem(my_cat_name)
    mmcif_cat.setValue(my_cat_val)
    return


def uniq_val(my_dict):
    """
    Makes the values of dictionary unique if k, v=[]
    """
    uni_dict = {}
    for key in my_dict:
        val = list(set(my_dict[key]))
        uni_dict[key] = sorted(val)
    return uni_dict


def transpose_list(my_data):
    np_array = np.array(my_data)
    transpose = np_array.T
    transpose_list = transpose.tolist()
    return transpose_list


def get_ent_chains(a, b, c, d):
    """
    Get a entity_chain dictionary, k=entity_id v= [asym_ids]
    """
    my_chain = {}
    res_list = {}
    for i in range(0, len(a)):
        my_chain.setdefault(int(a[i]), []).append(b[i])
        res_list.setdefault(int(a[i]), {}).setdefault(b[i], {}).setdefault(
            int(c[i]), {}
        )["mon_id"] = d[i]

    return uniq_val(my_chain), res_list


def make_opt(my_str):
    opt_str = my_str
    if isinstance(my_str, str) and my_str.strip() == "":
        opt_str = None

    return opt_str


def get_obs(my_chain, seq_id, obs_val):
    """Get observed residue from mmcif"""
    my_obs = {}
    for ch, i, j in zip(my_chain, seq_id, obs_val):
        my_obs.setdefault(ch, {})[int(i)] = j
    return my_obs


def get_mh_id(seq_id, mon_id, hetero):
    """
    Get mh_id - index for microheterogenity from mmcif
    """
    my_mh_id = {}
    boo = 1
    for i, j, k in zip(seq_id, mon_id, hetero):
        my_mh_id.setdefault(int(i), {})[j] = boo
        if k == "n":
            boo = 1
        else:
            boo = boo + 1

    return my_mh_id


def modify_atomsite(atom_site: dict[str, list], sifts_data):
    """
    Add db_name/db_num/db_res/db_acc in atomsite
    """

    for i in range(len(atom_site["label_seq_id"])):
        resnum = int(atom_site["label_seq_id"][i])
        asym = atom_site["label_asym_id"][i]
        ent = int(atom_site["label_entity_id"][i])

        db_name, unp_res, unp_acc, unp_num = None, None, None, None

        if (
            ent in sifts_data
            and asym in sifts_data[ent]
            and resnum in sifts_data[ent][asym]
        ):
            if len(sifts_data[ent][asym][resnum]) == 1:
                tar = sifts_data[ent][asym][resnum][0]

                db_name = "UNP"
                unp_res, unp_acc, unp_num = tar[2], tar[3], tar[4]

            else:
                tar1 = [
                    f"{item[2]}_{item[3]}_{item[4]}"
                    for item in sifts_data[ent][asym][resnum]
                ]
                tar1 = list(set(tar1))
                if len(tar1) == 1:
                    db_name = "UNP"
                    unp_res, unp_acc, unp_num = tar1[0].split("_")
                else:
                    logger.debug(
                        f"Chromophorore residue: {sifts_data[ent][asym][resnum]}"
                    )
                    # these are chromophores - where one pdb_seq_id
                    # residue may be mapped to more than one unp_seq_id
                    # these are ignored at the moments
                    pass

        atom_site.setdefault("pdbx_sifts_xref_db_name", []).append(db_name)
        atom_site.setdefault("pdbx_sifts_xref_db_acc", []).append(unp_acc)
        atom_site.setdefault("pdbx_sifts_xref_db_num", []).append(unp_num)
        atom_site.setdefault("pdbx_sifts_xref_db_res", []).append(unp_res)

    return atom_site


def expand_xref_seg_to_resi(my_res, xref_data):
    """expand xref_db segments mapping to every residue"""
    my_xref_db = ["Pfam", "CATH", "SCOP2"]
    exp_resi = {}
    for entity in sorted(my_res):
        for chain in my_res[entity]:
            # get boundaries for pfam,scop2,cath
            all_res = my_res[entity][chain].keys()
            if entity in xref_data and chain in xref_data[entity]:
                for xref_db in my_xref_db:
                    if xref_db in xref_data[entity][chain]:
                        for start in xref_data[entity][chain][xref_db]:
                            for region in xref_data[entity][chain][xref_db][start]:
                                end, acc, name, seg_id, inst_id = (
                                    region[0],
                                    region[1],
                                    region[2],
                                    region[3],
                                    region[4],
                                )
                                # print(start,end,acc,name)
                                xref_res = [
                                    item for item in all_res if start <= item <= end
                                ]
                                for resi in xref_res:
                                    exp_resi.setdefault(entity, {}).setdefault(
                                        chain, {}
                                    ).setdefault(resi, {}).setdefault(xref_db, {})[
                                        acc
                                    ] = [
                                        name,
                                        seg_id,
                                        inst_id,
                                    ]
    return exp_resi


def expand_unp_seg_to_resi(my_res, sifts_seg_inst):
    """expand unp segments mapping to every residue"""
    my_unp = {}
    for entity in sorted(my_res):
        for chain in my_res[entity]:
            # get boundaries for pfam,scop2,cath
            all_res = my_res[entity][chain].keys()
            if entity in sifts_seg_inst and chain in sifts_seg_inst[entity]:
                for start in sifts_seg_inst[entity][chain]:
                    for region in sifts_seg_inst[entity][chain][start]:
                        end, acc, seg_id, inst_id = (
                            region[0],
                            region[1],
                            region[2],
                            region[3],
                        )

                        my_unp_res = [item for item in all_res if start <= item <= end]
                        for resi in my_unp_res:
                            my_unp.setdefault(entity, {}).setdefault(
                                chain, {}
                            ).setdefault(resi, {})[acc] = [seg_id, inst_id]

    return my_unp


def get_xref_db(
    my_res, sifts_data, xref_data, sifts_seg_inst, my_obs, my_mh_id, my_mon_info
):
    """
    Merge sifts_res_csv and xref_res_data
    @data sorted here based on _poly_seq_schema
    """

    mega_list = [
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    ]

    expa_xref = expand_xref_seg_to_resi(my_res, xref_data)
    expa_sifts = expand_unp_seg_to_resi(my_res, sifts_seg_inst)
    # order in which xref_db are written
    my_xref_db = ["Pfam", "CATH", "SCOP2"]

    for entity in sorted(my_res):
        for chain in my_res[entity]:
            for res in sorted(my_res[entity][chain]):
                mon_id = my_res[entity][chain][res]["mon_id"]
                unp_res, unp_num, unp_acc = None, None, None
                res_type, mh_id, observed = None, None, None
                unp_seg_id, unp_inst_id = None, None
                db_acc, db_name, db, db_segment, db_instance = (
                    None,
                    None,
                    None,
                    None,
                    None,
                )
                observed = my_obs[chain][res]
                mh_id = my_mh_id[res][mon_id]
                boo = 1  # seq_id_instance/ordinal

                if (
                    entity in sifts_data
                    and chain in sifts_data[entity]
                    and res in sifts_data[entity][chain]
                ):
                    for tag in sifts_data[entity][chain][res]:
                        # mon_id        =   tag[0]
                        mon_oneLetter = make_opt(tag[1])
                        unp_res = make_opt(tag[2])
                        unp_acc = make_opt(tag[3])
                        unp_num = make_opt(tag[4])
                        res_type = make_opt(tag[5])
                        mh_id = make_opt(tag[6])
                        unp_seg_id, unp_inst_id = expa_sifts[entity][chain][res][
                            unp_acc
                        ]
                        if (
                            (
                                entity in expa_xref
                                and chain in expa_xref[entity]
                                and res not in expa_xref[entity][chain]
                            )
                            or (expa_xref == {})
                            or (entity not in expa_xref)
                            or (entity in expa_xref and chain not in expa_xref[entity])
                        ):
                            xx = [
                                entity,
                                chain,
                                boo,
                                res,
                                mon_id,
                                mon_oneLetter,
                                unp_res,
                                unp_num,
                                unp_acc,
                                unp_seg_id,
                                unp_inst_id,
                                res_type,
                                observed,
                                mh_id,
                                db_name,
                                db_acc,
                                db,
                                db_segment,
                                db_instance,
                            ]
                            [
                                my_list.append(my_val)
                                for my_list, my_val in zip(mega_list, xx)
                            ]
                            boo = boo + 1

                        else:
                            for xref_db in my_xref_db:
                                if (
                                    entity in expa_xref
                                    and chain in expa_xref[entity]
                                    and res in expa_xref[entity][chain]
                                    and xref_db in expa_xref[entity][chain][res]
                                ):
                                    # print(pfam[res])
                                    for acc in expa_xref[entity][chain][res][xref_db]:
                                        name, seg_id, inst_id = expa_xref[entity][
                                            chain
                                        ][res][xref_db][acc]
                                        xx = [
                                            entity,
                                            chain,
                                            boo,
                                            res,
                                            mon_id,
                                            mon_oneLetter,
                                            unp_res,
                                            unp_num,
                                            unp_acc,
                                            unp_seg_id,
                                            unp_inst_id,
                                            res_type,
                                            observed,
                                            mh_id,
                                            xref_db,
                                            acc,
                                            name,
                                            seg_id,
                                            inst_id,
                                        ]
                                        [
                                            my_list.append(my_val)
                                            for my_list, my_val in zip(mega_list, xx)
                                        ]
                                        boo = boo + 1

                else:
                    try:
                        mon_oneLetter = my_mon_info[entity][chain][res][mon_id]
                    except KeyError:
                        mon_oneLetter = None
                    for xref_db in my_xref_db:
                        if (
                            entity in expa_xref
                            and chain in expa_xref[entity]
                            and res in expa_xref[entity][chain]
                            and xref_db in expa_xref[entity][chain][res]
                        ):
                            # print(pfam[res])
                            for acc in expa_xref[entity][chain][res][xref_db]:
                                name, seg_id, inst_id = expa_xref[entity][chain][res][
                                    xref_db
                                ][acc]
                                xx = [
                                    entity,
                                    chain,
                                    boo,
                                    res,
                                    mon_id,
                                    mon_oneLetter,
                                    unp_res,
                                    unp_num,
                                    unp_acc,
                                    unp_seg_id,
                                    unp_inst_id,
                                    res_type,
                                    observed,
                                    mh_id,
                                    xref_db,
                                    acc,
                                    name,
                                    seg_id,
                                    inst_id,
                                ]
                                [
                                    my_list.append(my_val)
                                    for my_list, my_val in zip(mega_list, xx)
                                ]
                                boo = boo + 1

    return mega_list


def check_output_mmcif(file1, file2):
    """Checks if the file1 contents are still present in file2"""
    block_a = cif.read(str(file1)).sole_block()
    block_b = cif.read(str(file2)).sole_block()

    for cat in block_a.get_mmcif_category_names():
        a_cat_table = block_a.get_mmcif_category(cat)
        b_cat_table = block_b.get_mmcif_category(cat)
        if "_atom_site" in cat:
            common = set(a_cat_table.keys()).intersection(set(b_cat_table.keys()))
            for key in common:
                assert a_cat_table[key] == b_cat_table[key], f"{key} for {cat} differs"
        else:
            assert a_cat_table and b_cat_table and a_cat_table == b_cat_table


def check_sifts_mmcif(file1, file2, category_list):
    """Checks if the sifts category contents in file1 are still present in file2"""

    # infile = mmcifIO.CifFileReader(input="data", preserve_order=True)

    block_a = cif.read(str(file1)).sole_block()
    block_b = cif.read(str(file2)).sole_block()

    for cat in category_list:
        logger.debug("Category: {cat}")
        a_cat_table = block_a.get_mmcif_category(cat)
        b_cat_table = block_b.get_mmcif_category(cat)

        if "_atom_site" in cat:
            common = set(a_cat_table.keys()).intersection(set(b_cat_table.keys()))
            for key in common:
                assert a_cat_table[key] == b_cat_table[key]
        else:
            assert a_cat_table and b_cat_table and a_cat_table == b_cat_table
