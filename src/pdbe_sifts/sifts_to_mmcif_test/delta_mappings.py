import csv
import gzip
import os

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import get_next_release_date
from pdbe_sifts.sifts_to_mmcif.def_mmcif_cat import MAPPED_DB


def read_sifts_segments(block):
    """
    Read the segments information from sifts mmcif block
    and stores in a list and dictionary for further comparisons
    """
    data = {}
    rows = []
    cat = "_pdbx_sifts_unp_segments"
    mcat = block.get_mmcif_category(cat)
    if mcat:
        my_chain, my_unp_acc, my_unp_start, my_unp_end, my_pdb_start, my_pdb_end = (
            mcat["asym_id"],
            mcat["unp_acc"],
            mcat["unp_start"],
            mcat["unp_end"],
            mcat["seq_id_start"],
            mcat["seq_id_end"],
        )
        for chain, unp_acc, unp_start, unp_end, pdb_start, pdb_end in zip(
            my_chain, my_unp_acc, my_unp_start, my_unp_end, my_pdb_start, my_pdb_end
        ):
            rows.append(
                f"{chain}_UNP_{unp_acc}_{pdb_start}_{pdb_end}_{unp_start}_{unp_end}"
            )
            data.setdefault(chain, {}).setdefault("UNP", {}).setdefault(
                unp_acc, []
            ).append(f"{pdb_start}_{pdb_end}_{unp_start}_{unp_end}")
    cat = "_pdbx_sifts_xref_db_segments"
    mcat = block.get_mmcif_category(cat)
    if mcat:
        my_chain, my_db, my_db_acc, my_pdb_start, my_pdb_end = (
            mcat["asym_id"],
            mcat["xref_db"],
            mcat["xref_db_acc"],
            mcat["seq_id_start"],
            mcat["seq_id_end"],
        )
        for chain, db, db_acc, pdb_start, pdb_end in zip(
            my_chain, my_db, my_db_acc, my_pdb_start, my_pdb_end
        ):
            rows.append(f"{chain}_{db}_{db_acc}_{pdb_start}_{pdb_end}")
            data.setdefault(chain, {}).setdefault(db, {}).setdefault(db_acc, []).append(
                f"{pdb_start}_{pdb_end}__"
            )

    rows = sorted(rows)
    return data, rows


def get_delta_csv_suffix(csv_dir):
    """
    Program keeps changes from last 10(+1) releases.
    This function gets the suffix for each release and
    remove the oldest delta csv file  if more than 10.
    """
    my_cutoff = 10
    suff = str(get_next_release_date()).split("-", 1)[1]
    # file all the files with delta mappings
    all_file = sorted(
        [
            os.path.join(csv_dir, f)
            for f in os.listdir(csv_dir)
            if f.endswith(".csv.gz")
        ],
        key=os.path.getctime,
    )
    if len(all_file) >= my_cutoff:
        logger.info(f"All files: {all_file}")
        # delete the files which are extra
        to_del_file = all_file[: len(all_file) - my_cutoff]
        # print(to_del_file)
        for item in to_del_file:
            logger.info(f"Removing old delta_mappings {item}")
            os.remove(item)

    return suff


def write_delta_csv(del_map, sifts_delta_csv):
    """
    Write the delta mapping csv file
    """
    with gzip.open(sifts_delta_csv, "wt") as csvfile:
        del_file = csv.writer(
            csvfile,
            delimiter=",",
            quotechar="'",
            quoting=csv.QUOTE_ALL,
            lineterminator="\n",
        )
        for item in del_map:
            del_file.writerow(item)


class FindMappingChanges:
    """
    Returns DELTA_TYPE i.e.what has changed,
    delta_type can has three values-
        1.new_mapping
        2.obsoleted_mapping
        3.changed_segments
    Value returs is a list of x where x is
        x= [release_date,delta_type,pdb,chain,database,db_acc,old_segments,new_segments]
            where:
                database can be MAPPED_DB (from def_mmcif_cat)
                    current db = ['UNP','Pfam','CATH','SCOP2']

                old_segments/new_segments is None when not applicable
                    for eg new_mapping/obsoleted_mapping
                old_segments/new_segments is ";" seperated segments,
                where each segment is f{pdb_start}_{pdb_end}_{unp_start}_{unp_end}

    """

    def __init__(self, pdb, old_data, new_data, sifts_delta_csv):
        self.pdb = pdb
        self.old_data = old_data
        self.new_data = new_data
        self.sifts_delta_csv = sifts_delta_csv
        self.delta_lol = []
        self.rel_date = str(get_next_release_date())
        self.my_db = MAPPED_DB
        self.all_chains = sorted(self.old_data.keys() | self.new_data.keys())
        self.changed_db = []

    def old_vs_new(self, source_db, source_chain):
        """
        Compare the old and new data for a given database.
        It categories changes into three categories:
        1. obsoleted_mapping when the old_acc is not in present in new runs
        2. new_mappings when the acc is not seen in previous run
        3. changed_segments when the acc is same but mapped regions are different
        """
        old_val, new_val = [], []
        if source_db in self.old_data[source_chain]:
            old_val = self.old_data[source_chain][source_db].keys()
        if source_db in self.new_data[source_chain]:
            new_val = self.new_data[source_chain][source_db].keys()

        missing = [item for item in old_val if item not in new_val]
        new = [item for item in new_val if item not in old_val]
        comm = [item for item in new_val if item in old_val]
        # print(chain_id,old_val,new_val
        if missing:
            # print("obselted_mapping",source_chain,source_db,missing)
            self.handle_missing_acc(source_chain, source_db, missing)
        if new:
            self.handle_new_acc(source_chain, source_db, new)
        if comm:
            self.handle_common_acc(source_chain, source_db, comm)

        return

    def handle_missing_acc(self, source_chain, source_db, missing):
        """
        Marks the accession which are obsolted now
        """

        # print("obselted_mapping",source_chain,source_db,missing)
        for acc in missing:
            for seg in self.old_data[source_chain][source_db][acc]:
                # print("obsoleted_mapping",chain_id,source_db,acc,seg,"")
                self.changed_db.append(source_db)
                self.delta_lol.append(
                    [
                        self.rel_date,
                        "obsoleted_mapping",
                        self.pdb,
                        source_chain,
                        source_db,
                        acc,
                        seg,
                        None,
                    ]
                )

    def handle_new_acc(self, source_chain, source_db, new):
        """
        Marks the accession which are newly mapped
        """
        for acc in new:
            for seg in self.new_data[source_chain][source_db][acc]:
                # print("new_mapping",chain_id,source_db,acc,"",seg)
                self.changed_db.append(source_db)
                self.delta_lol.append(
                    [
                        self.rel_date,
                        "new_mapping",
                        self.pdb,
                        source_chain,
                        source_db,
                        acc,
                        None,
                        seg,
                    ]
                )

    def handle_common_acc(self, source_chain, source_db, comm):
        """
        Marks the accessions where segments or mapped regions have changed
        """
        for acc in comm:
            if (
                self.new_data[source_chain][source_db][acc]
                != self.old_data[source_chain][source_db][acc]
            ):
                old_seg = ";".join(self.old_data[source_chain][source_db][acc])
                new_seg = ";".join(self.new_data[source_chain][source_db][acc])
                # print("changed_segments",chain_id,source_db,acc,old_seg,new_seg)
                self.changed_db.append(source_db)
                self.delta_lol.append(
                    [
                        self.rel_date,
                        "changed_segments",
                        self.pdb,
                        source_chain,
                        source_db,
                        acc,
                        old_seg,
                        new_seg,
                    ]
                )

    def handle_lost_chain_mappings(self, source_chain):
        """
        Add details for cases where-
        mapping is obsoleted and complete chain is removed.
        """
        for db in self.old_data[source_chain]:
            for acc in self.old_data[source_chain][db]:
                seg = ";".join(self.old_data[source_chain][db][acc])
                # print("obselete_mapping",source_chain,db,acc,seg,"")
                self.changed_db.append(db)
                self.delta_lol.append(
                    [
                        self.rel_date,
                        "obsoleted_mapping",
                        self.pdb,
                        source_chain,
                        db,
                        acc,
                        seg,
                        None,
                    ]
                )

    def handle_new_chain_mappings(self, source_chain):
        """
        Add details for cases where-
        entire new chain is mapped.
        """

        for db in self.new_data[source_chain]:
            for acc in self.new_data[source_chain][db]:
                seg = ";".join(self.new_data[source_chain][db][acc])
                # print("new_mapping",source_chain,db,acc,"",seg)
                self.changed_db.append(db)
                self.delta_lol.append(
                    [
                        self.rel_date,
                        "new_mapping",
                        self.pdb,
                        source_chain,
                        db,
                        acc,
                        None,
                        seg,
                    ]
                )

    def find_mapping_changes(self):
        """
        Find all the changes in mapping and
        writes in csv file
        header for this file is:
        release_date,delta_type,pdb,chain,database,db_acc,old_segments,new_segments
        where old/new_segments are ";" seperated segments where,
        each segment is seperate by "_",{pdb_start}_{pdb_end}_{unp_start}_{unp_end}
        """

        for chain in self.all_chains:
            if chain in self.old_data:
                if chain in self.new_data:
                    for db in self.my_db:
                        self.old_vs_new(db, chain)
                else:
                    self.handle_lost_chain_mappings(chain)
            else:
                self.handle_new_chain_mappings(chain)

        if self.delta_lol:
            self.changed_db = sorted(set(self.changed_db))
            logger.info(f"Mapping changed for {self.pdb} {','.join(self.changed_db)}")
            write_delta_csv(self.delta_lol, self.sifts_delta_csv)
        return
