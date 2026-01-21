#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
import tempfile
import duckdb
from pathlib import Path

from gemmi import cif

from pdbe_sifts.base import pdbe_path
from pdbe_sifts.base.exceptions import EntryFailedException
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.parser import parse_with_base_parser
from pdbe_sifts.config import Config
from pdbe_sifts.sifts_to_mmcif import comm_utils, fetch_xref_data
from pdbe_sifts.sifts_to_mmcif.def_mmcif_cat import (
    NEW_MMCIF_CAT,
    PRIMARY_KEYS,
    SIFTS_ATOMSITE_ITEM,
    SIFTS_NEW_CAT,
)
from pdbe_sifts.sifts_to_mmcif.delta_mappings import (
    FindMappingChanges,
    get_delta_csv_suffix,
    read_sifts_segments,
)
from pdbe_sifts.sifts_to_mmcif.read_sifts_csv import get_unp_segments, get_unpres_mapping

conf = Config()


class NoSegmentsError(Exception):
    pass


class ExportSIFTSTommCIF():
    """
    Write the updated mmcif file with sifts_data
    """

    def __init__(
        self,
        output_path,
        cif_dir,
        sifts_csv_dir,
        duckdb_file,
        track_changes,
        prev_run_dir,
        delta_file,
    ):
        self.output_path = output_path
        self.cif_dir = cif_dir
        self.duckdb_file = duckdb_file
        self.sifts_csv_base = sifts_csv_dir
        self.conn = duckdb.connect(self.duckdb_file)

        self.track_changes = track_changes
        self.prev_run_dir = prev_run_dir
        self.entry_file_path = conf.lists.entries_all

        self.delta_file = delta_file


    def _check_mmcif_output(self, infile, outfile):
        logger.debug("Checking if the output file is more/equal to input file")

        if outfile.exists():
            input_size = infile.stat().st_size
            output_size = outfile.stat().st_size
            if output_size < input_size:
                logger.warning(
                    "Check output file- it's smaller than input) for "
                    f"{outfile} opsize vs insize --> {output_size} vs {input_size}"
                )
                comm_utils.check_output_mmcif(infile, outfile)

        else:
            raise EntryFailedException(f"Output {outfile} file not found for {outfile}")

    def _check_clean_mmcif(self, outfile):
        """
        Checking if the output clean mmcif file doesn't have null values in their primary keys
        """
        block = cif.read(str(outfile)).sole_block()
        category_list = [f"_{item}" for item in SIFTS_NEW_CAT]
        required_fields = PRIMARY_KEYS
        required_fields.setdefault("_atom_site", []).append(
            [
                "pdbx_sifts_xref_db_name",
                "pdbx_sifts_xref_db_acc",
                "pdbx_sifts_xref_db_num",
                "pdbx_sifts_xref_db_res",
            ]
        )
        for cat in category_list:
            cols = required_fields[cat]
            for col in cols:
                if any(cif.is_null(x) for x in block.find_values(f"{cat}.{col}")):
                    raise Exception(
                        f"Data error invalid value (? or .) for required item: {cat}.{col}"
                    )

    def _check_sifts_mmcif(self, sifts_only_cif: Path):
        logger.debug("Checking sifts mmcif output is not empty")
        if sifts_only_cif.exists():
            if sifts_only_cif.stat().st_size <= 40:
                raise EntryFailedException(f"Empty sifts output file {sifts_only_cif}")
        else:
            raise EntryFailedException(
                f"file not found sifts mmcif file {sifts_only_cif}"
            )

    def _process_cif_file(
        self, d_block, sifts_segments, sifts_seg_inst, sifts_res_csv, pdb_id
    ):
        # parsing the input mmcif file
        self.UpdateFlag = False
        self.UnpDataFlag = False

        sifts_res_info = {}
        # for each data block in mmcif file

        """
        Adding new category #1: pdbx_sifts_unp_segments
        @values coming here from sifts_segment_csv
        """

        self.UpdateFlag = self.write_data(
            d_block, "_pdbx_sifts_unp_segments", sifts_segments
        )
        # print(d_block.get_mmcif_category("_pdbx_sifts_unp_segments"))
        """
        Adding new category #2: pdbx_sifts_xref_db_segments
        @values coming here from database
        """
        # Get total entity_id,asymid for pdbentry, for sorting your data
        cat = d_block.get_mmcif_category("_pdbx_poly_seq_scheme")
        if cat:
            my_ent = cat["entity_id"]
            my_ch = cat["asym_id"]
            my_seq_id = cat["seq_id"]
            my_mon_id = cat["mon_id"]
            # get observed residues from mmcif _pdbx_poly_seq_scheme
            my_observed = cat["auth_seq_num"]
            my_observed = ["y" if item is not None else "n" for item in my_observed]
            my_obs_res = comm_utils.get_obs(my_ch, my_seq_id, my_observed)
            # get mh_id flag from mmcif
            hetero = cat["hetero"]
            my_mh_id = comm_utils.get_mh_id(my_seq_id, my_mon_id, hetero)

            all_ent, all_res = comm_utils.get_ent_chains(
                my_ent, my_ch, my_seq_id, my_mon_id
            )
            # xref_segments, xref_res_info = fetch_xref_data.get_xref_info(
            #     all_ent, pdb_id, self.conn
            # )
            xref_segments, xref_res_info = [], {}

            self.UpdateFlag = self.write_data(
                d_block, "_pdbx_sifts_xref_db_segments", xref_segments
            )

            """
            Adding new category #3: pdbx_sifts_xref_db
            showed best mapped uniprot accession
            @values coming here from sifts_residue_csv
            """
            if sifts_res_csv:
                res_cursor = None
            else:
                res_cursor = self.conn
            sifts_res_info, mon_id_info = get_unpres_mapping(
                pdb_id, sifts_res_csv, res_cursor
            )
            xref_resmap = comm_utils.get_xref_db(
                all_res,
                sifts_res_info,
                xref_res_info,
                sifts_seg_inst,
                my_obs_res,
                my_mh_id,
                mon_id_info,
            )
            xref_resmap = [item for item in xref_resmap if item != []]

            self.UpdateFlag = self.write_data(
                d_block, "_pdbx_sifts_xref_db", xref_resmap
            )

        else:
            logger.warning(f"The entry {pdb_id} is not a polymer")

        """
        Modifying atom site category #4:
        add db_name, db_accession, db_num, db_res
        modified only if sifts unp mapping present,else skipped
        """
        if sifts_res_info:
            self.UnpDataFlag = True
            self.UpdateFlag = True
            atom_site = d_block.get_mmcif_category("_atom_site")

            # Rename pdbe_label_seq_id
            atom_site = {
                "pdbx_label_index" if k == "pdbe_label_seq_id" else k: v
                for k, v in atom_site.items()
            }

            # add new items to _atom_site category
            atom_site = comm_utils.modify_atomsite(atom_site, sifts_res_info)

            d_block.set_mmcif_category("_atom_site", atom_site)
        else:
            if not self.UnpDataFlag:
                atom_site = d_block.get_mmcif_category("_atom_site")
                # Rename pdbe_label_seq_id
                atom_site = {
                    "pdbx_label_index" if k == "pdbe_label_seq_id" else k: v
                    for k, v in atom_site.items()
                }
                d_block.set_mmcif_category("_atom_site", atom_site)

        return self.UpdateFlag, d_block, self.UnpDataFlag

    def write_data(self, d_block, category, category_data):
        if category_data and len(NEW_MMCIF_CAT[category]) == len(category_data):
            data = dict(zip(NEW_MMCIF_CAT[category], category_data))
            d_block.set_mmcif_category(category, data)
            return True
        return False

    def delta_changes(self):
        # track delta mappings before writing the new file
        # reading sifts segments from current run
        self.sifts_delta_csv = None
        if self.track_changes and self.UnpDataFlag:
            base_sifts_dir = Path(self.sifts_updated_cif).parent
            delta_suff = get_delta_csv_suffix(base_sifts_dir)
            self.sifts_delta_csv = Path(
                Path(self.sifts_updated_cif).parent,
                f"{self.entry_id}_delta_{delta_suff}.csv.gz",
            )
            new_seg_data, new_seg_list = read_sifts_segments(self.d_block)
            # reading sifts segment from old/previous run
            # read from prev run dir given
            if self.prev_run_dir:
                old_sifts_csv_dir = pdbe_path.get_entry_dir(
                    self.entry_id, base_dir=self.prev_run_dir, suffix="sifts"
                )
                old_sifts_cif = Path(
                    Path(old_sifts_csv_dir), f"{self.entry_id}_sifts_only.cif.gz"
                )
            # look in the output folder as prev run dir not given
            else:
                old_sifts_cif = self.sifts_only_cif

            old_seg_data, old_seg_list = {}, []
            if old_sifts_cif.exists():
                old_block = cif.read(str(old_sifts_cif)).sole_block()
                old_seg_data, old_seg_list = read_sifts_segments(old_block)
            else:
                logger.warning(
                    f"Track changes failed for {self.entry_id}.\
                    Previous run file not found!Either a new entry or previous run failed!"
                )
            if old_seg_list == new_seg_list:
                logger.info(
                    f"No mappings changes tracked for {self.entry_id}.\
                    Skipping writing new mmcif files"
                )
                self.sifts_delta_csv = None
                return
            else:
                logger.info(f"Tracked mappings changed for {self.entry_id}")
                FindMappingChanges(
                    self.entry_id, old_seg_data, new_seg_data, self.sifts_delta_csv
                ).find_mapping_changes()

    def process_entry(self, entry_id):
        """
        @requiremnts clean_cif_file,seg_csv,res_csv file
        @database connection
        """
        sifts_csv_dir = None
        if self.sifts_csv_base:
            sifts_csv_dir = pdbe_path.get_entry_dir(
                entry_id, base_dir=self.sifts_csv_base, suffix="sifts"
            )
        self.sifts_updated_cif = Path(
            pdbe_path.get_sifts_updated_cif(entry_id, base_dir=self.output_path)
        )
        self.sifts_only_cif = Path(
            Path(self.sifts_updated_cif).parent, f"{entry_id}_sifts_only.cif.gz"
        )
        Path(self.sifts_updated_cif).parent.mkdir(parents=True, exist_ok=True)

        if sifts_csv_dir:
            seg_cursor = None
            sifts_seg_csv = Path(sifts_csv_dir, f"{entry_id}_seg.csv.gz")
            if not sifts_seg_csv.exists():
                logger.error(f"{sifts_seg_csv} does not exist")

            sifts_res_csv = Path(sifts_csv_dir, f"{entry_id}_res.csv.gz")
            if not sifts_seg_csv.exists():
                logger.error(f"{sifts_res_csv} does not exist")
        else:
            sifts_seg_csv, sifts_res_csv = None, None
            seg_cursor = self.conn

        sifts_segments, sifts_seg_inst = get_unp_segments(
            entry_id, sifts_seg_csv, seg_cursor
        )
        sifts_segments = [item for item in sifts_segments if item != []]

        # uses pdb_updated.cif.gz
        orig_updated_cif = Path(pdbe_path.get_updated_cif(entry_id, self.cif_dir))

        if not Path(orig_updated_cif).exists():
            raise FileNotFoundError(orig_updated_cif)

        self.d_block = cif.read(str(orig_updated_cif)).sole_block()

        self.UpdateFlag, self.d_block, self.UnpDataFlag = self._process_cif_file(
            self.d_block, sifts_segments, sifts_seg_inst, sifts_res_csv, entry_id
        )

        # track delta mappings before writing the new file
        # reading sifts segments from current run
        self.entry_id = entry_id
        self.delta_changes()

        # writing new mmcif file
        doc = cif.Document()
        doc.add_copied_block(self.d_block)
        with gzip.open(self.sifts_updated_cif, "wt") as f:
            f.write(doc.as_string(cif.Style.Pdbx))

        # writing a snipped file :
        # removes the original categories keeping only sifts categories
        # for next-gen-archive, we need snipped file always for pdbx_label_index

        doc = cif.Document()
        block = doc.add_new_block(entry_id.upper())
        # copy required sifts categories
        cat_to_copy = SIFTS_NEW_CAT

        for cat in cat_to_copy:
            data = self.d_block.get_mmcif_category(f"_{cat}")
            if not data:
                if cat == "_pdbx_sifts_unp_segments":
                    raise NoSegmentsError("No segments written. This shouldn't happen")
                logger.warning(f"Category _{cat} not written out")
                continue
            block.set_mmcif_category(f"_{cat}.", data)

        # take only limited items from atom_site
        atom_cat = block.get_mmcif_category("_atom_site")
        if self.UnpDataFlag:
            for key in list(atom_cat.keys()):
                if key not in SIFTS_ATOMSITE_ITEM:
                    del atom_cat[key]
        else:
            # copy the atom_site id and pdbx_label_index
            for key in list(atom_cat.keys()):
                if key not in ["id", "pdbx_label_index"]:
                    del atom_cat[key]
        block.set_mmcif_category("_atom_site.", atom_cat)

        # Write to temp file. This is needed to avoid
        # corrupting the original file in case of errors
        # while writing the file
        with tempfile.TemporaryDirectory() as d:
            temp_name = Path(d, f"{entry_id}_sifts_only.cif.gz")
            with gzip.open(temp_name, "wt") as f:
                f.write(doc.as_string(cif.Style.Pdbx))
                f.flush()
            # Copy temp file to sifts_only_cif
            shutil.copyfile(temp_name, self.sifts_only_cif)

        self._check_sifts_mmcif(self.sifts_only_cif)
        del self.d_block
        del block
        del doc

        if self.UpdateFlag:
            logger.info(f"SIFTS mapping found for {self.entry_id}")

        else:
            logger.info(f"No mappings in SIFTS for {self.entry_id}")
        # need to delete the object to read file again

        self._check_mmcif_output(orig_updated_cif, self.sifts_updated_cif)

        self._check_clean_mmcif(self.sifts_updated_cif)
        logger.info(
            f"""
            Input CIF: {orig_updated_cif}
            Output files:
            UPDATED_MMCIF: {self.sifts_updated_cif}
            SIFTS_ONLY: {self.sifts_only_cif}
            DELTA_MAPPINGS: {self.sifts_delta_csv}
        """
        )
        return (orig_updated_cif, self.sifts_updated_cif, self.sifts_only_cif)

    def after_job_end(self):
        """Close the database connection."""
        close_db_connection(self.conn)

    def post_run(self):
        # works only in batch mode
        # generate a list of pdb entries and what db mapping changed
        if self.track_changes and self.delta_file and self.mode == Modes.BATCH:
            logger.info("Batching all pdb entries whose mapping changed from the logs")
            cmd = f"""
            grep Mapping {self.log_dir}/* |cut -d"]" -f2 |cut -d" " -f6,7 >{self.delta_file}"""
            subprocess.check_call(cmd, shell=True)
            logger.info(f"DELTA_MAPPING_LIST: {self.delta_file}")


def run():
    parser = argparse.ArgumentParser("Given a clean mmcif file, add sifts_data")
    parser.add_argument(
        "-o",
        "--output-dir",
        help="directory where output mmcif files will be stored",
        default=conf.sifts_to_mmcif.clean_output_path,
        # action=OrcAction,
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input-cif-dir",
        # action=OrcAction,
        required=True,
        default=conf.location.work.data_entry_dir,
        help="input cif base directory",
    )

    parser.add_argument(
        "-s",
        "--sifts-csv-dir",
        # action=OrcAction,
        required=False,
        help="sifts_xref_residue csv file, if not given data taken from database",
    )

    parser.add_argument(
        "-d",
        "--db-conn-str",
        required=True,
        # action=OrcAction,
        # envvar="ORC_PDBEREAD_DB_CONN",
        help="Database connection string",
    )

    parser.add_argument(
        "-T",
        "--no-track-changes",
        required=False,
        action="store_true",
        default=False,
        help="""
        Compares the mapping from previous run and track changes.
        Default is True. Store the delta_mappings_csv from last 10 releases
        """,
    )

    parser.add_argument(
        "-p",
        "--prev-run-dir",
        required=False,
        # action=OrcAction,
        help="""
        Compare sifts_only.mmcif from here.
        If not given, it looks in the output directory for previous run file,
        it expects data in the format it writes data i.e. given_dir/suff_pdb/pdb/sifts/.
        If no file is found, it assumes it is a new csv file
        """,
    )

    parser.add_argument(
        "-l",
        "--delta-sifts-file",
        required=False,
        # action=OrcAction,
        default=conf.lists.sifts_mapping_changes,
        help="""
        Directory where delta_sifts.list is generated (works only in batch mode).
        This file has pdb entries which has there mapping changed last week.
        It contains pdb, comma seperated list of db for which mapping changed.
        """,
    )

    args = parse_with_base_parser(parser)

    track_changes = not args.no_track_changes
    c = ExportSIFTSTommCIF(
        args.output_dir,
        args.input_cif_dir,
        args.sifts_csv_dir,
        args.db_conn_str,
        track_changes,
        args.prev_run_dir,
        args.delta_sifts_file,
    )
    c.process_entry(args.entry)


if __name__ == "__main__":
    run()
