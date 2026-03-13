import itertools
from pathlib import Path

from gemmi import cif

from pdbe_sifts.config import load_config

conf = load_config()

# mmCIF categories used in SIFTS
CATEGORIES = [
    "_pdbx_poly_seq_scheme",
    "_entity_poly",
    "_pdbx_struct_mod_residue",
    "_struct_ref_seq_dif",
    "_struct_ref",
    "_struct_ref_seq",
    "_entity_src_nat",
    "_entity_src_gen",
    "_pdbx_entity_src_syn",
    "_entity",
    "_pdbx_database_status",
    "_pdbx_audit_revision_history",
]

def extract_table(block, search_list):
    """
    Produces a Gemmi table based on a list of parsed column names in the _atom_site.
    loop.

    Input: Gemmi block, list of search terms
    Returned: Gemmi table

    :param mmcif: Contents of (updated) mmCIF file
    :type mmcif: gemmi.cif.Block
    :param search_list: Loop terms to extract
    :type search_list: list[str]
    :return: Values for searched loop terms
    :rtype: gemmi.cif.Table
    """

    table = block.find(
        "_atom_site.",
        search_list,
    )

    return table


class NotAPolyPeptide(Exception):
    pass


class mmCIF:
    """Docstring for mmCIF."""

    def __init__(self, pdbid, chem_comp_dict, cif_file: str):
        self.pdbid = pdbid
        self.fname = cif_file
        if not Path(self.fname).exists():
            raise FileNotFoundError(f"The mmcif file {self.fname} does not exists.")

        block = cif.read(self.fname).sole_block()
        self.poly_seq = block.get_mmcif_category("_pdbx_poly_seq_scheme")

        if not self.poly_seq:
            raise NotAPolyPeptide(
                "it is not a polypeptide! no poly_seq, only hetatms presents"
            )
        self.cc = chem_comp_dict
        self.entity_poly = block.get_mmcif_category("_entity_poly")
        self.mod_residue = block.get_mmcif_category("_pdbx_struct_mod_residue")
        self.seq_dif = block.get_mmcif_category("_struct_ref_seq_dif")
        # self.ref_seq = block.get_mmcif_category("struct_ref_seq")
        self.struct_ref = block.get_mmcif_category("_struct_ref")
        self.struct_ref_seq = block.get_mmcif_category("_struct_ref_seq")
        self.src_nat = block.get_mmcif_category("_entity_src_nat")
        self.src_gen = block.get_mmcif_category("_entity_src_gen")
        self.src_syn = block.get_mmcif_category("_pdbx_entity_src_syn")
        self.entity = block.get_mmcif_category("_entity")
        self.pdb_status = block.get_mmcif_category("_pdbx_database_status")
        self.rev_date = block.get_mmcif_category("_pdbx_audit_revision_history")
        self.features = self.__get_features()
        del block

    def get_unp(self, chain):
        """Fetches uniprot accessions from CIF for a particular entity, where available."""
        unps: list[str] = []

        if self.struct_ref:
            entity = self.get_entity_id(chain)
            entities = self.struct_ref["entity_id"]
            names = self.struct_ref["db_name"]
            accs = self.struct_ref["pdbx_db_accession"]

            for ent, acc, name in zip(entities, accs, names):
                if name == "UNP" and ent == entity:
                    unps.append(acc)

        return unps

    def get_ranges(self, chain, acc):
        if self.struct_ref_seq is None:
            return None

        chains = self.struct_ref_seq["pdbx_strand_id"]
        accs = self.struct_ref_seq["pdbx_db_accession"]
        begins = self.struct_ref_seq["db_align_beg"]
        ends = self.struct_ref_seq["db_align_end"]

        out = []

        if isinstance(chains, list):
            for idx, c in enumerate(chains):
                if c == chain and accs[idx] == acc:
                    out.append([int(begins[idx]), int(ends[idx])])
        else:
            if chains == chain and accs == acc:
                out = [[int(begins), int(ends)]]

        return out

    def __get_tax_helper(self, cat, item_name, entity):
        out = None

        if cat:
            entities = cat["entity_id"]

            if not isinstance(entities, list):
                entities = [entities]
                is_list = False
            else:
                is_list = True

            for idx, ent in enumerate(entities):
                if entity == ent:
                    if is_list:
                        out = cat[item_name][idx]
                        break
                    else:
                        out = cat[item_name]
                        break

        if out is None:
            return None

        # Comma-separated list
        if "," in out:
            out = out.split(",")[0]

        try:
            return int(out)
        except ValueError:
            return None

    def get_tax(self, chain):
        entity = self.get_entity_id(chain)

        tax = self.__get_tax_helper(self.src_nat, "pdbx_ncbi_taxonomy_id", entity)

        if tax is not None:
            return tax

        tax = self.__get_tax_helper(
            self.src_gen, "pdbx_gene_src_ncbi_taxonomy_id", entity
        )

        if tax is not None:
            return tax

        tax = self.__get_tax_helper(self.src_syn, "ncbi_taxonomy_id", entity)

        return tax

    def get_chains(self):
        """TODO: Docstring for get_chains.

        @param f TODO
        @return: TODO

        """

        chains = self.poly_seq["pdb_strand_id"]
        out = {}

        for idx, c in enumerate(chains):
            if c in out:
                continue
            out[c] = (
                self.poly_seq["entity_id"][idx],
                self.poly_seq["asym_id"][idx],
            )

        return out

    def is_poly(self, entity_id):
        # Check if the entity is a polypeptide
        if isinstance(entity_id, int):
            entity_id = str(entity_id)

        cat = self.entity_poly

        for idx, entity in enumerate(cat["entity_id"]):
            if entity == entity_id:
                if isinstance(cat["type"], list):
                    return "peptide" in cat["type"][idx].lower()
                else:
                    return "peptide" in cat["type"].lower()

        return False

    def get_pdb_status(self):
        return self.pdb_status["status_code"][0]

    def get_rev_date(self):
        return self.rev_date["revision_date"]

    def get_type(self, chain, n_residue):
        if not self.seq_dif:
            return None

        n_residue = str(n_residue)

        chains = self.seq_dif["pdbx_pdb_strand_id"]

        if not isinstance(chains, list):
            chains = [chains]
            is_list = False
        else:
            is_list = True

        if chain not in chains:
            return None

        for idx, c in enumerate(chains):
            if c == chain:
                if is_list:
                    if n_residue == self.seq_dif["seq_num"][idx]:
                        try:
                            return self.seq_dif["details"][idx].capitalize()
                        except AttributeError:
                            return None
                else:
                    if n_residue == self.seq_dif["seq_num"]:
                        try:
                            return self.seq_dif["details"].capitalize()
                        except AttributeError:
                            return None

        return None

    def get_original_residue(self, chain, n_residue):
        if self.seq_dif is None:
            return None

        n_residue = str(n_residue)

        chains = self.seq_dif["pdbx_pdb_strand_id"]

        if chain not in chains:
            return None

        for idx, c in enumerate(chains):
            if c == chain:
                if n_residue == self.seq_dif["seq_num"][idx]:
                    threeL = self.seq_dif["db_mon_id"][idx]
                    return (threeL, self.cc.get(threeL))

        return None

    def get_sequence(self, entity_id):
        # TODO: is this overwritting sequence positions by going through all
        # the chains in the entity?  that is not a problem in terms of accuracy
        # but it is not optimal in terms of computational time
        # mon = poly_seq["mon_id"]
        if isinstance(entity_id, int):
            entity_id = str(entity_id)

        entities = self.poly_seq["entity_id"]
        seq = {}

        for idx, e in enumerate(entities):
            if e != entity_id:
                continue

            seq_id = int(self.poly_seq["seq_id"][idx])

            if seq_id not in list(seq.keys()):
                threeL = self.poly_seq["mon_id"][idx]

                seq[seq_id] = self.cc.get(threeL)

        return "".join([seq[key] for key in sorted(seq.keys())])

    def __get_parent_id(self, three):
        label = self.mod_residue["label_comp_id"]
        parent = self.mod_residue["parent_comp_id"]

        for lbl, prnt in zip(label, parent):
            if lbl == three:
                return prnt

        return None

    def __get_features(self):
        if not self.seq_dif:
            return {}

        features = []
        out = {}

        chains = self.seq_dif["pdbx_pdb_strand_id"]
        seq = self.seq_dif["seq_num"]
        details = self.seq_dif["details"]

        features = [
            (c, int(s), d)
            for c, s, d in zip(chains, seq, details)
            if s not in ("?", ".", None, False)
        ]

        s = sorted(features, key=lambda x: (x[0], x[1], x[2]))

        for key, groups in itertools.groupby(s, key=lambda x: (x[0], x[2])):
            if key[0] not in out:
                out[key[0]] = []

            idx = 0
            groups = list(groups)

            start = groups[0][1]
            prev = start
            idx = 1

            while idx < len(groups):
                if groups[idx][1] != prev + 1:
                    out[key[0]].append((key[1], (start, prev)))
                    start = groups[idx][1]

                prev = groups[idx][1]
                idx += 1

            out[key[0]].append((key[1], (start, prev)))

        return out

    def get_entity_id(self, chain):
        chains = self.poly_seq["pdb_strand_id"]

        for idx, c in enumerate(chains):
            if c == chain:
                return self.poly_seq["entity_id"][idx]

    def get_residues(self, chain):
        chains = self.poly_seq["pdb_strand_id"]

        residues = {}

        for idx, c in enumerate(chains):
            if c != chain:
                continue

            n = int(self.poly_seq["seq_id"][idx])

            observed = True if self.poly_seq["auth_seq_num"][idx] else False
            auth_n = self.poly_seq["pdb_seq_num"][idx] if observed else None
            auth_ins = self.poly_seq["pdb_ins_code"][idx] if observed else None

            threeL = self.poly_seq["mon_id"][idx]
            oneL = self.cc.get(threeL)
            mh = self.poly_seq["hetero"][idx] == "y"
            rtype = self.get_type(c, n)
            oneL_original = None
            threeL_original = None

            if rtype in ("Engineered mutation", "Conflict"):
                out = self.get_original_residue(c, n)
                if out is not None:
                    threeL_original, oneL_original = out

            # auth_n, auth_ins, oneL, threeL, rtype, observed, MH
            tmp = (
                auth_n,
                auth_ins,
                oneL,
                threeL,
                rtype,
                observed,
                oneL_original,
                threeL_original,
            )

            if mh:
                if n in residues:
                    residues[n].append(tmp + (len(residues[n]) + 1,))
                else:
                    residues[n] = [tmp + (1,)]
            else:
                residues[n] = tmp + (1,)

        return sorted(residues.items())

    def get_ec(self, entity_id):
        # Get the EC for the entity
        if isinstance(entity_id, int):
            entity_id = str(entity_id)

        cat = self.entity
        out = None

        if cat.get("pdbx_ec") is None:
            return []

        for idx, entity in enumerate(cat["id"]):
            if entity == entity_id:
                if isinstance(cat["pdbx_ec"], list):
                    out = cat["pdbx_ec"][idx]
                else:
                    out = cat["pdbx_ec"]
                break

        if out is not None:
            return [x.strip() for x in out.split(",") if x not in ("?", ".")]

        return []

    def get_ters_coordinates_of_residue(self, chain_id, res_num, ters_lab_num):
        # function from Joseph Ellaway
        cter_lab = ters_lab_num["C"][0]
        nter_lab = ters_lab_num["N"][0]
        cter_num = int(ters_lab_num["C"][1])
        nter_num = int(ters_lab_num["N"][1])
        block = cif.read(self.fname).sole_block()
        res_num = str(res_num)

        mmcif_table = extract_table(
            block,
            [
                "group_PDB",
                "id",
                "type_symbol",
                "label_atom_id",
                "label_asym_id",
                "label_seq_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "pdbx_PDB_model_num",
            ],
        )

        if len(mmcif_table) == 0:
            # Could not parse mmCIF file
            logger.error(
                f"Gemmi block {block} for chain {chain_id} does not contain valid "
                "pdbx_sifts_xref_db_num (UniProt sequence ID) column"
            )

            raise TypeError(
                "Updated mmCIF file is missing SIFTs column headers needed for clustering. "
                "Please repeat clustering step for this UniProt segment once mmCIFs have "
                "been updated."
            )

        # Loop over Table and make checks
        first_model = mmcif_table[0][9]
        previous_residue = mmcif_table[0][5]
        atom_ordinal = 1
        return_dict = {"N": [], "C": []}
        for row in mmcif_table:
            if row[5] != previous_residue:
                previous_residue = row[5]
                atom_ordinal = 1
            if row[9] != first_model:
                break
            atom_lab = row[3]
            # Check: Atom is C or N in protein residue, with occupancy and in correct chain
            if (
                row[4] == chain_id
                and row[5] == res_num
                and (atom_lab == cter_lab or atom_lab == nter_lab)
                and (atom_ordinal == cter_num or atom_ordinal == nter_num)
            ):
                # Add Cartesian x,y,z
                return_dict[row[2]] = [
                    float(row[6]),
                    float(row[7]),
                    float(row[8]),
                ]
            atom_ordinal += 1
        return return_dict

