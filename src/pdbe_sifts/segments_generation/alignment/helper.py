from enum import Enum
from multiprocessing.dummy import Pool
from typing import NamedTuple

import tqdm
from funcy.debug import log_durations

from pdbe_sifts.base.log import logger
from pdbe_sifts.segments_generation import alignment
from pdbe_sifts.mmcif.entry import Entry
# from segments_gen.uniref90_pkl import NF90Coverage, NF90TaxID
from pdbe_sifts.unp.unp import UNP
from pdbe_sifts.unp.unp import get_unp_object

NF_COVERAGE = 0.7
N_PROC = 64
STEP_SIZE = 2000


class SplitAccessionError(Exception):
    pass


class SMapping(NamedTuple):
    accession: str
    range_start: int
    range_stop: int


class SiftsMode(Enum):
    ISOFORM = 1
    NF90 = 2


def get_accession(entry: Entry, acc: str) -> UNP:
    if acc in entry.accessions:
        return entry.accessions[acc]

    unp = get_unp_object(acc)
    if unp:
        entry.accessions[unp.accession] = unp

    return unp


def check_range(real_ranges, pdb, unp):
    for p, u in zip(pdb, unp):
        if (p, u) not in real_ranges:
            return False
    return True


def fmt_ranges(ranges):
    out = []

    for r in ranges:
        if isinstance(r, str):
            out.append([int(x) for x in r.split("-")])
        else:
            out.append(r)

    return out


def overlapping(ranges):
    ranges = fmt_ranges(ranges)

    for idx, r1 in enumerate(ranges):
        for r2 in ranges[idx + 1 :]:
            if r2[0] <= r1[0] <= r2[1] or r2[0] <= r1[1] <= r2[1]:
                return True

    return False


def out_of_order(ranges):
    ranges = fmt_ranges(ranges)

    for idx, r1 in enumerate(ranges):
        for r2 in ranges[idx + 1 :]:
            if r2[0] <= r1[0]:
                return True

    return False


def unroll_map(m):
    chain = m[0]

    if len(m) > 1:
        accs = eval(m[1])
        ranges = eval(m[2])
    else:
        accs = None
        ranges = None

    return chain, accs, ranges


# Currently:
#   GAP_OPEN = 10
#   GAP_EXTEND = 2
#
# Therefore, gap won't open if (AT LEAST!):
# G > E / 2 - 5
#
# G = gap length
# E = extra aligned sequence
#
# This method identifies manually annotated large gaps so we can treat them
# as overlapping mappings (forcing extra alignments, one per range)
def large_gap(ranges):
    if len(ranges) < 2:
        return False

    # If ranges come from mmCIF: array of ints
    # If ranges come from AD_DBREF (redis): array of strings
    if isinstance(ranges[0], str):
        iranges = [list(map(int, x.split("-"))) for x in ranges]
    else:
        iranges = ranges

    gap = iranges[1][0] - iranges[0][1] - 1
    extra = iranges[1][1] - iranges[1][0] + 1

    return gap > 3 * (extra / 2 - 5)


class EntryMapping:
    def __init__(
        self,
        entry: Entry,
        chain,
        chain_mapping: list[SMapping],
        nf90_mode=False,
        NFT=None,
        NFC=None,
    ):
        self.entry = entry
        self.chain = chain
        self.chain_mapping = chain_mapping
        self.nf90_mode = nf90_mode
        self.chain_obj = self.entry.chains[self.chain]
        self.repeated_acc = False
        self.accs: list[str] = []
        self.ranges: list[tuple[int, int]] = []
        self.NFT: NF90TaxID = NFT
        self.NFC: NF90Coverage = NFC

    def _get_accessions(self):
        # Validate the provided mappings
        for mapping in self.chain_mapping:
            unp = get_accession(self.entry, mapping.accession)
            if unp:
                logger.info(
                    f"Valid mapping: {unp.accession} {mapping.range_start}-{mapping.range_stop}"
                )
                self.accs.append(unp.accession)
                self.ranges.append((mapping.range_start, mapping.range_stop))

        # If no valid mappings were provided, try to get them from the mmCIF
        if not self.accs:
            cif_accs = self.entry.mmcif.get_unp(self.chain)
            for acc in cif_accs:
                unp = get_accession(self.entry, acc)
                if unp:
                    self.accs.append(unp.accession)

            if self.accs:
                logger.info(f"Accessions come from mmCIF: {self.accs}")

        # If we have mappings but no ranges, try to get them from the mmCIF
        if self.accs and not self.ranges:
            self.ranges = self.entry.mmcif.get_ranges(self.chain, self.accs[0])
            logger.info(f"Ranges come from mmCIF: {self.ranges}")

    def set_chain_accessions(self):
        self._get_accessions()
        if not self.accs:
            logger.warning("The mmCIF doesnt have a valid UniProt accession")
            return False

        logger.debug(f"Accessions [{self.chain}]: {self.accs}")
        try:
            self.chain_obj = self.entry.chains[self.chain]
        except KeyError:
            logger.warning(
                f"The chain {self.chain} was not found in the mmCIF  or it is not a polypeptide"
            )
            return False

        if self.nf90_mode:
            if not self.check_unp_coverage():
                return False
        return True

    @log_durations(logger.debug)
    def process(self):
        self.seq_pdb = self.chain_obj.alignment_sequence
        self.check_repeated_accession()

        logger.info(f"Processing {self.chain}: {self.accs}")

        for acc in set(self.accs):
            self.acc = acc
            isoforms = {}
            unp = get_accession(self.entry, acc)
            if not unp:
                continue
            self.unp = unp
            # isoform, score, length
            self.chain_obj.best[unp.accession] = (None, 0.0, 0)

            # Add the canonical
            self.chain_obj.canonicals.append(unp.accession)

            # Don't store isoforms if the chain is a chimera
            if self.chain_obj.is_chimera:
                isoforms[unp.accession] = unp.sequence
            elif self.nf90_mode:
                isoforms = self.get_nf90_isoforms(unp, isoforms)
            else:
                for iso, seq in list(unp.getAllIsoforms().items()):
                    isoforms[iso] = seq
                    self.entry.accessions[iso] = unp

            # create a process pool that uses all cpus
            with Pool(N_PROC) as pool:
                # call the function for each item in parallel, get results as tasks complete
                my_list = list(isoforms.items())
                list(
                    tqdm.tqdm(
                        pool.imap_unordered(
                            self.process_each_isoform, my_list, chunksize=STEP_SIZE
                        )
                    )
                )

        self.chain_obj.generate_residue_maps()
        logger.info("Segments generated")

    def process_each_isoform(self, row):
        iso, seq = row
        # Don't repeat the canonical
        if (
            not self.nf90_mode
            and seq == self.unp.sequence
            and iso != self.unp.accession
        ):
            logger.info(f"{iso} is the canonical. Skipping...")
            return iso

        all_alns = alignment.do_alignment_lalign36(seq, self.seq_pdb)
        alns = []
        # We want as many alignments as times the accession is repeated.
        # Alignments sorted from best to worst
        if self.repeated_acc:
            for _ in range(self.accs.count(self.acc)):
                try:
                    alns.append(next(all_alns))
                except StopIteration:
                    break
        else:
            try:
                alns.append(next(all_alns))
            # the alignment might be impossible
            # (e.g. an isoform which doesn't have the relevant sequence fragment)
            except StopIteration:
                pass
        self.process_alignments(self.unp, iso, alns)
        return iso

    def process_alignments(self, unp, iso, alns):
        # Loop through as many alignments as accessions in the list
        # if it is only one accession which repeats
        for al in alns:
            pdb_start, pdb_end = (al[1]._al_start, al[1]._al_stop)
            unp_start, unp_end = (al[0]._al_start, al[0]._al_stop)

            self.chain_obj.mod_after_alignment(al)

            logger.debug(alignment.annotate_alignment(al[0].seq, al[1].seq))

            logger.debug("PDB: [%d-%d]" % (pdb_start, pdb_end))
            logger.debug("UNP: [%d-%d]" % (unp_start, unp_end))

            pdb_ranges = [(pdb_start, pdb_end)]
            unp_ranges = [(unp_start, unp_end)]

            self.chain_obj.mappings.setdefault(iso, []).append(
                (pdb_ranges, unp_ranges, al)
            )
            identity = alignment.get_identity(al[0]._seq, al[1]._seq)
            score = alignment.get_score(al[0]._seq, al[1]._seq)

            self.chain_obj.scores[iso] = (identity, score)

            self.chain_obj.seg_scores.setdefault(iso, {})[str(pdb_ranges)] = identity

            self._update_best_mapping(unp, iso, al, score)

    def get_nf90_isoforms(self, unp, isoforms):
        logger.debug(f"fetching NF90 isoforms for {unp}")
        for acc in self.NFT.get_nf90(unp.accession):
            unp_nf90 = get_accession(self.entry, acc)

            if not unp_nf90:
                continue

            acc = unp_nf90.accession if "-" not in acc else acc

            if acc not in unp_nf90.seq_isoforms:
                unp_nf90.seq_isoforms[acc] = unp_nf90.sequence

            if acc not in self.entry.accessions:
                self.entry.accessions[acc] = unp_nf90

            isoforms[acc] = unp_nf90.sequence
        return isoforms

    def get_valid_accession(self, acc):
        unp = get_accession(self.entry, acc)
        if not unp:
            accs = self.entry.mmcif.get_unp(self.chain)
            if not accs:
                logger.warning(
                    f"The mmCIF doesnt have a UniProt accession: {self.entry.pdbid}"
                )
                self.chain_obj.is_chimera = False
                self.repeated_acc = False
                return

            if len(accs) > 1:
                logger.error(
                    "It is a chimera and at least one of the accessions is not valid anymore."
                )
                self.chain_obj.is_chimera = False
                self.repeated_acc = False
                return

            logger.warning(
                f"The accession {acc} is not valid. Got {accs[0]} from mmCIF"
            )
            if acc == accs[0]:
                logger.warning(
                    f"CIF has same invalid obsolete accession {acc}. Doing nothing"
                )
                return

            unp = get_accession(self.entry, accs[0])

        return unp

    def _update_best_mapping(self, unp, iso, al, score):
        try:
            ad_dbref_acc = self.chain_mapping[0].accession
        except Exception:
            logger.warning(
                f"Could not get the AD_DBREF accession for {self.entry.pdbid} {self.chain}"
            )
            ad_dbref_acc = None
        best_score = self.chain_obj.best[unp.accession][1]
        best_seq_len = self.chain_obj.best[unp.accession][2]
        if (
            score > best_score
            or (score == best_score and len(al[1]._seq) > best_seq_len)
            or (score == best_score and len(al[1]._seq) == best_seq_len)
            and (iso in self.chain_obj.canonicals or iso == ad_dbref_acc)
        ):
            self.chain_obj.best[unp.accession] = (
                iso,
                score,
                len(al[1]._seq),
            )

    def check_unp_coverage(self):
        logger.debug("Checking coverage")
        unp = get_accession(self.entry, self.accs[0])

        if not unp:
            return False

        if self.repeated_acc or self.chain_obj.is_chimera:
            logger.warning("It is a chimera so we skip it for UniRef90")
            return False
        coverage = self.NFC.get_coverage(self.entry.pdbid, self.chain, unp.accession)
        if coverage < NF_COVERAGE:
            logger.warning(
                f"Coverage {coverage} < {NF_COVERAGE} so we skip it for UniRef90"
            )
            return False
        return True

    def check_repeated_accession(self):
        if len(self.accs) > 1:
            # One accession that repeats
            if len(set(self.accs)) == 1:
                if out_of_order(self.ranges) or overlapping(self.ranges):
                    self.repeated_acc = True
                    logger.warning("Same accession repeats!")
                # Large gap go hereƒ
                elif large_gap(self.ranges):
                    # Calculate gap and extra aligned sequence and,
                    # if large enough, treat as a repeated_acc
                    self.repeated_acc = True
                    logger.warning("Large gap detected")
            else:
                logger.info(f"Chain {self.chain} is Chimera!")
                self.chain_obj.is_chimera = True
