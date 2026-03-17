#!/usr/bin/env python3
"""
Module to manage entities and entries from an mmCIF file.

An entity represents a unique sequence in a PDB structure.
Multiple chains can share the same entity (e.g. homodimers).
"""

from itertools import groupby
from multiprocessing.dummy import Pool
from operator import itemgetter

import tqdm
from Bio.Seq import Seq

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import get_cpu_count

from . import mmcif_helper
from .chain import Chain
from .residue import Residue

N_PROC = get_cpu_count()
STEP_SIZE = 2000


class Entity:
    """
    Represents a biological entity in an mmCIF file.

    An entity corresponds to a unique sequence (_entity_poly in mmCIF).
    Multiple chains (auth_asym_id) can belong to the same entity
    (for example, in a homodimer, chains A and B share the same entity).

    Attributes:
        mmcif (mmcif_helper.mmCIF): Parsed mmCIF object
        pdbid (str): PDB identifier (e.g. "1ABC")
        entity_id (str): Entity identifier (e.g. "1", "2", etc.)
        sequence (str): Primary sequence in one-letter code
        auth_asym_ids (List[str]): List of chain identifiers (e.g. ["A", "B"])
        struct_asym_ids (List[str]): List of internal structural identifiers
        ec (List[str]): Associated EC (Enzyme Commission) numbers
        tax_id (Optional[int]): NCBI taxonomy identifier
        parent (object): Reference to parent Entry object
        chains (Dict[str, Chain]): Dictionary of Chain objects for this entity
    """

    def __init__(
        self,
        mmcif: mmcif_helper.mmCIF,
        pdbid: str,
        entity_id: str,
        parent: object | None = None,
    ):
        """
        Initialize an entity from mmCIF data.

        Args:
            mmcif: mmCIF object containing parsed data
            pdbid: 4-character PDB identifier
            entity_id: Entity identifier
            sequence: Primary sequence in one-letter code
            parent: Reference to parent Entry object (optional)
        """
        self.mmcif = mmcif
        self.pdbid = pdbid
        self.entity_id = entity_id
        self.sequence = self.get_entity_sequence()
        self.parent = parent
        self.best = {}
        self.scores = {}
        self.seg_scores = {}
        self.canonicals = []
        self.is_chimera = False

        # alignment info
        # residues are the same to each chain belonging to the same entity
        # however in chain X they can be observed and not in chain Y
        self.residues: list[Residue] = []
        self.alignment_sequence: list[str] = []
        self.expression_tag_start = None

        # Lists of chains belonging to this entity
        self.auth_asym_ids: list[str] = []
        self.struct_asym_ids: list[str] = []

        # References to Chain objects for this entity
        self.chains: dict[str, Chain] = {}

        # Biological metadata
        self.ec: list[str] = []
        self.tax_id: int | None = None

        # mapping information
        self.mappings: dict = {}
        self.residue_maps = {}
        self.segments = {}

        # Flag indicating whether this entity should be skipped
        self.skip: bool = False

        # Initialization
        self._load_chains_info()
        self._load_metadata()
        self._load_residues()
        self._gen_alignment_sequence()

        logger.debug(
            f"Entity {entity_id} created for {pdbid} with {len(self.auth_asym_ids)} chain(s)"
        )

    def get_entity_sequence(self):
        return self.mmcif.get_sequence(self.entity_id)

    def _load_chains_info(self) -> None:
        """
        Load information about chains belonging to this entity.
        Note: Chain objects will be added later by Entry.
        """
        all_chains = self.mmcif.get_chains()
        # ==> {'auth_asym_id': ('entity', 'struct_asym')}

        for auth_asym_id, chain_data in all_chains.items():
            chain_entity_id = chain_data[0]
            struct_asym_id = chain_data[1]

            if chain_entity_id == self.entity_id:
                self.auth_asym_ids.append(auth_asym_id)
                self.struct_asym_ids.append(struct_asym_id)

        logger.debug(f"Entity {self.entity_id}: chains found = {self.auth_asym_ids}")

    def add_chain(self, auth_asym_id: str, chain: Chain) -> None:
        """
        Add a Chain object to this entity.

        Args:
            auth_asym_id: Chain identifier
            chain: Chain object
        """
        self.chains[auth_asym_id] = chain
        logger.debug(f"Chain {auth_asym_id} added to Entity {self.entity_id}")

    def _load_residues(self) -> None:
        """
        Load entity residues from mmCIF.

        Residues are extracted from a single representative chain
        (all chains of an entity share the same biological sequence).

        Observation (present / missing) is tracked per chain.
        """
        if not self.auth_asym_ids:
            logger.warning(f"Entity {self.entity_id}: no chains found")
            return

        # Storage: observed per chain, per seq_id
        self._observed_by_chain: dict[str, dict[int, bool]] = {}

        # --- 1. Load residues from the first chain (entity-level definition)
        ref_chain = self.auth_asym_ids[0]

        for n, v in self.mmcif.get_residues(ref_chain):
            if isinstance(v, list):
                v = v[0]

            residue = Residue(
                n,
                None,  # auth_n (chain level not entity)
                None,  # auth_ins (chain level not entity)
                v[2],  # oneL
                v[3],  # threeL
                v[4],  # rtype
                None,  # observed (chain level not entity)
                v[6],  # oneL_original
                v[7],  # threeL_original
                v[8],  # mh index
            )

            self.residues.append(residue)

        # --- 2. Collect observed flags for ALL chains
        for chain_id in self.auth_asym_ids:
            self._observed_by_chain[chain_id] = {}

            for n, v in self.mmcif.get_residues(chain_id):
                if isinstance(v, list):
                    v = v[0]
                self._observed_by_chain[chain_id][n] = v[5]

        logger.debug(
            f"Entity {self.entity_id}: loaded {len(self.residues)} residues "
            f"from chain {ref_chain}, observed states collected for "
            f"{len(self.auth_asym_ids)} chains"
        )

    def is_residue_observed(self, seq_id: int, chain_id: str) -> bool:
        """
        Return whether a residue is observed in a given chain.

        Args:
            seq_id: entity sequence position (1-based)
            chain_id: auth_asym_id

        Returns:
            True if observed, False otherwise
        """
        try:
            return self._observed_by_chain[chain_id].get(seq_id, False)
        except AttributeError:
            return False

    def _gen_alignment_sequence(self) -> None:
        """
        Generate the alignment sequence with modifications to simplify alignment.

        Applied modifications:
        - Insertions, Linkers, Expression tags → 'J'
        - Mutations, Conflicts, Cloning artifacts → original residue if available
        - Chromophores → expanded into their components
        """
        if self.alignment_sequence:
            return

        for r in self.residues:
            if isinstance(r, list):
                r = r[0]

            # Mark the start of the expression tag
            if r.rtype == "Expression tag" and self.expression_tag_start is None:
                self.expression_tag_start = r.n

            # Mask insertions, linkers, and expression tags
            if r.rtype in ("Insertion", "Linker", "Expression tag"):
                self.alignment_sequence.append("J")

            # Use original residue for mutations/conflicts
            elif r.rtype in ("Engineered mutation", "Conflict", "Cloning artifact"):
                self.alignment_sequence.append(r.oneL_original if r.oneL_original else "X")

            # Expand chromophores
            elif r.is_chromophore:
                self.alignment_sequence.extend(list(r.oneL))

            else:
                self.alignment_sequence.append(r.oneL)

        logger.debug(
            f"Entity {self.entity_id}: alignment sequence generated "
            f"({len(self.alignment_sequence)} positions)"
        )

    def mod_after_alignment(self, al):
        seq = al[1].seq
        out = ""

        i = 0

        while i < (al[1]._al_start - 1):
            r = self.residues[i]

            if isinstance(r, list):
                r = r[0]

            if len(r.oneL) > 1:
                al[1]._al_start -= len(r.oneL) - 1
                al[1]._al_stop -= len(r.oneL) - 1
                # i -= (len(r.oneL) - 1)

            i += 1

        i = al[1]._al_start - 1
        j = 0

        while j < len(seq):
            s = seq[j]

            if s == "-":
                out += s
                j += 1
                continue

            if isinstance(self.residues[i], list):
                r = self.residues[i][0]
            else:
                r = self.residues[i]

            out += r.oneL

            i += 1
            j += len(r.oneL)

        al[1].seq = Seq(out)
        return out

    # generate the residue_maps between the chain and each isoform
    def get_each_resmap(self, row):
        iso, mappings = row
        residue_map = {}
        # Process in reverse so the best mappings is
        # overriding the rest
        for m in mappings[::-1]:
            al = m[2]
            pdb_seq = al[1]._seq
            pdb_start = al[1]._al_start
            pdb_stop = al[1]._al_stop
            unp_seq = al[0]._seq
            unp_start = al[0]._al_start
            unp_stop = al[0]._al_stop
            pdb_i = 0
            unp_i = 0
            pdb_shift = 0
            unp_shift = 0
            while pdb_i + pdb_start <= pdb_stop and unp_i + unp_start <= unp_stop:
                pdb_r = pdb_seq[pdb_i + pdb_shift]
                unp_r = unp_seq[unp_i + unp_shift]
                r = self.residues[pdb_i + pdb_start - 1]
                if isinstance(r, list):
                    r = r[0]
                # Remove undesired mappings just in case they were included
                if (
                    unp_r != "-"
                    and pdb_r != "-"
                    and r.rtype in ("Insertion", "Linker", "Expression tag")
                ):
                    pdb_i += 1
                    unp_i += 1
                    continue
                if len(r.oneL) > 1 and pdb_r != "-" and unp_r != "-":
                    residue_map[pdb_i + pdb_start] = []
                    for x in range(len(r.oneL)):
                        residue_map[pdb_i + pdb_start].append(unp_i + unp_start + x)
                    unp_i += 1
                    pdb_i += 1
                elif pdb_r != "-" and unp_r != "-":
                    residue_map[pdb_i + pdb_start] = unp_i + unp_start
                    pdb_i += 1
                    unp_i += 1
                elif pdb_r == "-":
                    unp_i += 1
                    pdb_shift += 1
                elif unp_r == "-":
                    pdb_i += 1
                    unp_shift += 1
                if len(r.oneL) > 1:
                    pdb_shift += len(r.oneL) - 1
                    unp_i += len(r.oneL) - 1

        self.residue_maps[iso] = residue_map
        return iso, residue_map

    # generate the residue_map between the chain and each isoform
    def generate_residue_maps(self):
        # create a process pool that uses all cpus
        with Pool(N_PROC) as pool:
            # call the function for each item in parallel, get results as tasks complete
            my_list = list(self.mappings.items())
            list(tqdm.tqdm(pool.imap_unordered(self.get_each_resmap, my_list, chunksize=STEP_SIZE)))

        # remove the residues which map to more than one accession
        # keeping the ones that benefit continuity
        if self.is_chimera:
            self.__overlapping_residues()

        # generate the segments once we have the residue maps
        self.__segments_from_residues()

    # remove the residues which map to more than one accession
    # keeping the ones that benefit continuity
    # (only for chimeras)
    def __overlapping_residues(self):
        iso_list = list(self.residue_maps.keys())

        for idx, iso1 in enumerate(iso_list):
            for iso2 in iso_list[idx + 1 :]:
                maps1 = self.residue_maps[iso1]
                maps2 = self.residue_maps[iso2]

                for key in list(maps1.keys()):
                    if key in maps2:
                        # The mapping without conflict has preference
                        pdb_r = self.sequence[key - 1]
                        try:
                            unp1_r = self.parent.accessions[iso1].seq_isoforms[iso1][maps1[key] - 1]
                            unp2_r = self.parent.accessions[iso2].seq_isoforms[iso2][maps2[key] - 1]

                            # If one is a conflict and the other one is not
                            if unp1_r != unp2_r and pdb_r in (unp1_r, unp2_r):
                                if pdb_r == unp2_r:
                                    del maps1[key]
                                else:
                                    del maps2[key]

                                continue
                        except TypeError:
                            # maps1[key] and/or maps2[key] are lists
                            # (the residue is a chromophore)
                            pass

                        # Otherwise try to keep continuity
                        if key - 1 in list(maps1.keys()):
                            del maps2[key]
                        else:
                            del maps1[key]

    def __group_elements(self, lst):
        ranges = []

        for _, g in groupby(enumerate(lst), lambda i_x: i_x[0] - i_x[1][0]):
            group = list(map(itemgetter(1), g))
            ranges.append(((group[0][0], group[-1][0]), (group[0][1], group[-1][1])))

        # get only the min/max if there is a 1 to many mapping
        for idx, r in enumerate(ranges):
            if isinstance(r[1][0], list) or isinstance(r[1][1], list):
                start = min(r[1][0]) if isinstance(r[1][0], list) else r[1][0]
                end = max(r[1][1]) if isinstance(r[1][1], list) else r[1][1]

                ranges[idx] = (r[0], (start, end))

        return ranges

    def __overlapping(self, r1, r2):
        return r2[0] <= r1[0] <= r2[1] or r2[0] <= r1[1] <= r2[1]

    def __segments_from_residues(self):
        for iso, maps in list(self.residue_maps.items()):
            self.segments[iso] = self.__group_elements(
                (x, y) for x, y in sorted(maps.items(), key=itemgetter(0))
            )

        # remove overlapping segments for chimeras.
        # we create a final_set of non-overlapping segments and
        # mark all the ones overlapping with it for removal
        # by inserting them in to_remove
        if self.is_chimera:
            final_set = []
            to_remove = {}

            for iso, maps in list(self.segments.items()):
                for m, u in maps:
                    for r in final_set:
                        if self.__overlapping(m, r):
                            to_remove.setdefault(iso, []).append((m, u))
                            break

                    if m not in list(to_remove.values()):
                        final_set.append(m)

            for key, val in list(to_remove.items()):
                for seg in val:
                    self.segments[key].remove(seg)

    def _load_metadata(self) -> None:
        """
        Load biological metadata (EC numbers, taxonomy).
        """
        # EC numbers
        self.ec = self.mmcif.get_ec(self.entity_id)

        # Taxonomy – use the first chain
        if self.auth_asym_ids:
            self.tax_id = self.mmcif.get_tax(self.auth_asym_ids[0])

    def get_uniprot_accessions(self) -> list[str]:
        """
        Retrieve UniProt accessions for all chains of the entity.

        Returns:
            List of UniProt accessions (duplicates removed)
        """
        accessions = []
        for chain_id in self.auth_asym_ids:
            unps = self.mmcif.get_unp(chain_id)
            accessions.extend(unps)
        return list(set(accessions))

    def is_enzyme(self) -> bool:
        """
        Check whether the entity is an enzyme (has EC numbers).

        Returns:
            True if the entity has at least one EC number
        """
        return len(self.ec) > 0

    def __repr__(self) -> str:
        """Concise representation of the entity."""
        return (
            f"Entity(pdbid='{self.pdbid}', entity_id='{self.entity_id}', "
            f"chains={self.auth_asym_ids}, length={len(self.sequence)})"
        )

    def __str__(self) -> str:
        """Detailed representation of the entity."""
        lines = [
            f"Entity {self.entity_id} from PDB {self.pdbid}",
            f"Chains: {', '.join(self.auth_asym_ids)}",
            f"Sequence length: {len(self.sequence)} residues",
        ]

        if self.ec:
            lines.append(f"EC numbers: {', '.join(self.ec)}")

        if self.tax_id:
            lines.append(f"Taxonomy ID: {self.tax_id}")

        return "\n".join(lines)

    def __len__(self) -> int:
        """Return the sequence length."""
        return len(self.sequence)
