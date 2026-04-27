import itertools
from pathlib import Path

from gemmi import cif

from pdbe_sifts.base.exceptions import NotAPolyPeptide
from pdbe_sifts.base.log import logger

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


def extract_table(block: cif.Block, search_list: list[str]) -> cif.Table:
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


class mmCIF:
    """Parser and accessor for mmCIF data used by the SIFTS pipeline.

    Loads the mmCIF categories required for segment generation on
    construction and exposes typed accessor methods for chains, residues,
    UniProt cross-references, taxonomy IDs, and PDB metadata.
    """

    def __init__(self, pdbid: str, chem_comp_dict, cif_file: str) -> None:
        """Initialise by reading and caching all required mmCIF categories.

        Args:
            pdbid: PDB identifier (e.g. ``"1abc"``).
            chem_comp_dict: A :class:`ChemCompMapping` instance used to
                convert three-letter residue codes to one-letter codes.
            cif_file: Path to the ``.cif`` or ``.cif.gz`` file.

        Raises:
            FileNotFoundError: If *cif_file* does not exist.
            NotAPolyPeptide: If the file contains no
                ``_pdbx_poly_seq_scheme`` category.
        """
        self.pdbid = pdbid
        self.fname = cif_file
        if not Path(self.fname).exists():
            raise FileNotFoundError(
                f"The mmcif file {self.fname} does not exists."
            )

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

    def get_unp(self, chain: str) -> list[str]:
        """Fetches uniprot accessions from CIF for a particular entity, where available."""
        unps: list[str] = []

        if self.struct_ref:
            entity = self.get_entity_id(chain)
            entities = self.struct_ref["entity_id"]
            names = self.struct_ref["db_name"]
            accs = self.struct_ref["pdbx_db_accession"]

            for ent, acc, name in zip(entities, accs, names, strict=False):
                if name == "UNP" and ent == entity:
                    unps.append(acc)

        return unps

    def get_ranges(self, chain: str, acc: str) -> list[list[int]] | None:
        """Return the UniProt alignment ranges for a chain–accession pair.

        Reads ``_struct_ref_seq`` to find the database alignment begin/end
        positions for the given chain and UniProt accession.

        Args:
            chain: Author chain ID.
            acc: UniProt accession to query.

        Returns:
            List of ``[begin, end]`` integer pairs, one per aligned region.
            Returns ``None`` if ``_struct_ref_seq`` is absent, or an empty
            list if the chain/accession pair has no recorded ranges.
        """
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

    def __get_tax_helper(
        self, cat: dict | None, item_name: str, entity: str
    ) -> int | None:
        """Extract a taxonomy ID from a source category for a given entity.

        Args:
            cat: Parsed mmCIF category dict (e.g. ``self.src_nat``).
                May be ``None`` if the category is absent from the file.
            item_name: Name of the field that holds the taxonomy ID within
                *cat* (e.g. ``"pdbx_ncbi_taxonomy_id"``).
            entity: Entity ID string to match against ``cat["entity_id"]``.

        Returns:
            Integer taxonomy ID, or ``None`` if not found or not parseable.
        """
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

    def get_tax(self, chain: str) -> int | None:
        """Return the NCBI taxonomy ID for the source organism of a chain.

        Tries ``_entity_src_nat``, then ``_entity_src_gen``, then
        ``_pdbx_entity_src_syn`` in that order.

        Args:
            chain: Author chain ID.

        Returns:
            Integer NCBI taxonomy ID, or ``None`` if unavailable in all
            three source categories.
        """
        entity = self.get_entity_id(chain)

        tax = self.__get_tax_helper(
            self.src_nat, "pdbx_ncbi_taxonomy_id", entity
        )

        if tax is not None:
            return tax

        tax = self.__get_tax_helper(
            self.src_gen, "pdbx_gene_src_ncbi_taxonomy_id", entity
        )

        if tax is not None:
            return tax

        tax = self.__get_tax_helper(self.src_syn, "ncbi_taxonomy_id", entity)

        return tax

    def get_chains(self) -> dict[str, tuple]:
        """Return all polypeptide chains in the structure.

        Returns:
            Dict mapping author chain ID → ``(entity_id, asym_id)`` for the
            first occurrence of each unique chain in ``_pdbx_poly_seq_scheme``.
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

    def is_poly(self, entity_id: str | int) -> bool:
        """Return ``True`` if the entity is a polypeptide.

        Args:
            entity_id: Entity identifier (string or integer).

        Returns:
            ``True`` when ``_entity_poly.type`` contains ``"peptide"`` for
            the given entity, ``False`` otherwise.
        """
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

    def get_pdb_status(self) -> str:
        """Return the PDB release status code (e.g. ``"REL"``)."""
        return self.pdb_status["status_code"][0]

    def get_rev_date(self) -> list[str]:
        """Return all revision dates from ``_pdbx_audit_revision_history``."""
        return self.rev_date["revision_date"]

    def get_type(self, chain: str, n_residue: int) -> str | None:
        """Return the residue annotation type for a specific chain position.

        Looks up ``_struct_ref_seq_dif`` for a match on chain and residue
        number and returns the capitalised ``details`` field (e.g.
        ``"Engineered mutation"``, ``"Insertion"``, ``"Conflict"``).

        Args:
            chain: Author chain ID.
            n_residue: 1-based sequence position.

        Returns:
            Capitalised annotation string, or ``None`` if no annotation is
            found or the category is absent.
        """
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

    def get_original_residue(self, chain: str, n_residue: int) -> tuple | None:
        """Return the original residue before an engineered mutation or conflict.

        Args:
            chain: Author chain ID.
            n_residue: 1-based sequence position.

        Returns:
            ``(three_letter_code, one_letter_code)`` tuple of the original
            residue, or ``None`` if not found.
        """
        if self.seq_dif is None:
            return None

        n_residue = str(n_residue)

        chains = self.seq_dif["pdbx_pdb_strand_id"]

        if chain not in chains:
            return None

        for idx, c in enumerate(chains):
            if c == chain and n_residue == self.seq_dif["seq_num"][idx]:
                threeL = self.seq_dif["db_mon_id"][idx]
                return (threeL, self.cc.get(threeL))

        return None

    def get_sequence(self, entity_id: str | int) -> str:
        """Return the one-letter amino acid sequence for an entity.

        Iterates over ``_pdbx_poly_seq_scheme`` and collects the one-letter
        code for each unique sequence position belonging to *entity_id*.
        Positions are resolved via the chemical component dictionary.

        Args:
            entity_id: Entity identifier (string or integer).

        Returns:
            One-letter sequence string ordered by ascending sequence position.

        Note:
            When the entity spans multiple chains, all chains are visited but
            each sequence position is only recorded once (the first occurrence
            wins), so the result is the entity-level canonical sequence.
        """
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

    def __get_parent_id(self, three: str) -> str | None:
        """Look up the parent residue code for a modified residue.

        Searches ``_pdbx_struct_mod_residue`` for *three* in the
        ``label_comp_id`` column and returns the corresponding
        ``parent_comp_id``.

        Args:
            three: Three-letter code of the modified residue to look up
                (e.g. ``"MSE"``).

        Returns:
            Three-letter code of the parent standard residue (e.g. ``"MET"``),
            or ``None`` if *three* is not listed as a modified residue.
        """
        label = self.mod_residue["label_comp_id"]
        parent = self.mod_residue["parent_comp_id"]

        for lbl, prnt in zip(label, parent, strict=False):
            if lbl == three:
                return prnt

        return None

    def __get_features(self) -> dict:
        """Parse ``_struct_ref_seq_dif`` into a per-chain feature index.

        Groups consecutive residues that share the same chain and annotation
        type into contiguous ranges.  The result is stored in ``self.features``
        and is used downstream to identify insertions, mutations, and other
        sequence differences.

        Returns:
            Dict mapping author chain ID → list of ``(detail_str, (start, end))``
            tuples, where *start* and *end* are 1-based sequence positions.
            Returns an empty dict when ``_struct_ref_seq_dif`` is absent.
        """
        if not self.seq_dif:
            return {}

        features = []
        out = {}

        chains = self.seq_dif["pdbx_pdb_strand_id"]
        seq = self.seq_dif["seq_num"]
        details = self.seq_dif["details"]

        features = [
            (c, int(s), d)
            for c, s, d in zip(chains, seq, details, strict=False)
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

    def get_entity_id(self, chain: str) -> str | None:
        """Return the entity ID for a given author chain identifier.

        Searches ``_pdbx_poly_seq_scheme`` for the first row whose
        ``pdb_strand_id`` matches *chain* and returns the corresponding
        ``entity_id``.

        Args:
            chain: Author chain ID (e.g. ``"A"``).

        Returns:
            Entity ID string (e.g. ``"1"``), or ``None`` if *chain* is not
            found in the poly seq scheme.
        """
        chains = self.poly_seq["pdb_strand_id"]

        for idx, c in enumerate(chains):
            if c == chain:
                return self.poly_seq["entity_id"][idx]

    def get_residues(self, chain: str) -> list:
        """Return all residues for a chain as sorted (seq_id, data) pairs.

        Reads ``_pdbx_poly_seq_scheme`` and builds a dict keyed by integer
        sequence position.  Each value is a tuple of:
        ``(auth_n, auth_ins, oneL, threeL, rtype, observed, oneL_original,
        threeL_original, mh_index)``.  For microheterogeneity positions
        (``hetero == "y"``), the value is a list of such tuples (one per
        alternate conformer).

        Args:
            chain: Author chain ID whose residues should be returned.

        Returns:
            List of ``(seq_id, residue_data)`` tuples sorted by ascending
            sequence position.
        """
        chains = self.poly_seq["pdb_strand_id"]

        residues = {}

        for idx, c in enumerate(chains):
            if c != chain:
                continue

            n = int(self.poly_seq["seq_id"][idx])

            observed = bool(self.poly_seq["auth_seq_num"][idx])
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

    def get_ec(self, entity_id: str | int) -> list[str]:
        """Return the EC numbers associated with an entity.

        Reads the ``pdbx_ec`` field from the ``_entity`` category.  Handles
        comma-separated values and filters out placeholder characters
        (``"?"`` and ``"."``).

        Args:
            entity_id: Entity identifier (string or integer).

        Returns:
            List of EC number strings (e.g. ``["3.4.21.4"]``).  Returns an
            empty list when no EC data is present or the entity is not found.
        """
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

    def get_ters_coordinates_of_residue(
        self, chain_id: str, res_num: int, ters_lab_num: dict
    ) -> dict:
        """Return the Cartesian coordinates of the N- and C-terminal atoms of a residue.

        Reads ``_atom_site`` from the mmCIF file and locates the specific atom
        rows for the requested residue by matching chain, residue number, atom
        label, and atom ordinal.  Only the first model is considered.

        Args:
            chain_id: Internal ``label_asym_id`` (struct_asym) of the chain.
            res_num: ``label_seq_id`` (1-based entity sequence position) of the
                residue whose terminal-atom coordinates are requested.
            ters_lab_num: Dict with keys ``"C"`` and ``"N"``, each mapping to a
                ``(atom_label, atom_ordinal)`` pair identifying the C-terminal
                and N-terminal atoms respectively (e.g.
                ``{"C": ("C", "1"), "N": ("N", "1")}``).

        Returns:
            Dict with keys ``"N"`` and ``"C"``, each containing a
            ``[x, y, z]`` list of float coordinates.  A key's value is an
            empty list if the atom was not found.

        Raises:
            TypeError: If the ``_atom_site`` table is empty (mmCIF is missing
                the required columns).
        """
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
                and (atom_lab in (cter_lab, nter_lab))
                and (atom_ordinal in (cter_num, nter_num))
            ):
                # Add Cartesian x,y,z
                return_dict[row[2]] = [
                    float(row[6]),
                    float(row[7]),
                    float(row[8]),
                ]
            atom_ordinal += 1
        return return_dict
