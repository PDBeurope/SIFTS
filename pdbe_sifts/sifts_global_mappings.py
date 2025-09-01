#!/usr/bin/env python3

import argparse
import io
import csv
import json
from timeit import default_timer as timer
from pathlib import Path
import pickle

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.global_mappings.mmseqs_search import MmSearch
from pdbe_sifts.global_mappings.blastp import BlastP
from pdbe_sifts.global_mappings.global_mappings_parser import GlobMappingsParser
from pdbe_sifts.base.utils import get_date, make_path
from pdbe_sifts.global_mappings.scoring_function import get_ranked_mappings

class SiftsGlobalMappings():
    def __init__(
        self,
        cif_file,
        out_dir,
        db_file,
        tool = 'mmseqs',
        out_global_mappings = None,
        threads = 1,
    ):
        """Generate structure - sequence mappings.

        Aligns the sequences for each chain of the entry and aligns them against a pre-formated database.
        Also saved the mappings into the out_global_mappings file in the directory out_dir.

        Args:
            cif_file (path to file): path to the structure file (PDBx/mmCIF format required).
            out_dir (path to directory): location where the outputs will be saved.
            out_global_mappings (path to file): file name where to save global mappings results (default: out_dir/cif_file/global_mappings.csv)
            tool (str): blastp or mmseqs (default: mmseqs)
            db_file (path to file): the DB file to run against (in mmseqs format or blast format.)
        """

        self.cif_file = cif_file
        self.out_dir = Path(out_dir)
        self.fasta_files_path = self.out_dir / 'fasta_files'
        self.unp_dir = self.out_dir / 'unp_files'
        self.fasta_files_path.mkdir(parents=True, exist_ok=True)
        self.unp_dir.mkdir(parents=True, exist_ok=True)
        if not out_global_mappings:
            entry_name = Path(cif_file).stem
            standard_name = f'global_mappings_{entry_name}.csv'
            out_global_mappings = str(Path(out_dir) / entry_name / standard_name)
        self.out_global_mappings = out_global_mappings
        self.tool = tool
        self.db_file = db_file
        self.used_cif_categories = [
            "entity_poly",
            "pdbx_struct_mod_residue",
            "pdbx_poly_seq_scheme",
            "struct_ref_seq_dif",
            "struct_ref",
            "struct_ref_seq",
            "entity_src_nat",
            "entity_src_gen",
            "pdbx_entity_src_syn",
            "entity",
            "pdbx_database_status",
            "pdbx_audit_revision_history",
        ]
        self.entry = None
        self.result_file_path = {}
        self.mappings = {}
        self.ranked_mappings = {}
        self.threads = threads
    
    def generate_fasta(self, entity_seq_tax_dict, entry_name):
        merged_fasta = []
        for ent, seq_tax in entity_seq_tax_dict.items():
            seq = seq_tax[0]
            tax_id = seq_tax[1]
            header = f'>pdb|{entry_name}-{ent}|{tax_id}'
            content = f'{header}\n{seq}\n'
            merged_fasta.append(content)
        tmp_fasta_path = self.fasta_files_path / f'tmp_{entry_name}.fasta'
        with open(tmp_fasta_path, 'w') as f:
            f.write("".join(merged_fasta))
        return tmp_fasta_path

    def mmseqs_search(self, id, fake_fasta):
        now = get_date()
        output_path = make_path(self.out_dir, id, 'mmseqs', f'hits_{id}.tsv', now)
        tmp_fold = Path(output_path).parent / f'tmp_{id}'
        tmp_fold.mkdir(parents=True, exist_ok=True)
        result = MmSearch(fake_fasta, self.db_file, output_path, tmp_fold, self.threads)
        result.run()
        self.result_file_path[id] = output_path
    
    def blastp_search(self, id, fasta_path):
        output_path = make_path(self.out_dir, id, 'blastp', f'hits_{id}.json')
        blastp_search = BlastP(fasta_path, self.db_file, output_path, threads = self.threads)
        blastp_search.run()
        self.result_file_path[id] = output_path

    def search(self, id, tmp_fasta_path):
        match self.tool:
            case 'mmseqs':
                self.mmseqs_search(id, tmp_fasta_path)
            case 'blastp':
                self.blastp_search(id, tmp_fasta_path)
            case _:
                raise ValueError(f"Unsupported tool: {self.tool}")

    def process(self):
        logger.info("Processing [%s]" % self.cif_file)
        start_cif = timer()
        entry_name = Path(self.cif_file).stem
        self.entry = Entry(entry_name, self.cif_file)
        entity_seq_tax = self.entry.get_entity_seq_tax()
        tmp_fasta_path = self.generate_fasta(entity_seq_tax, entry_name)
        self.search(entry_name, tmp_fasta_path)
        self.mappings = GlobMappingsParser(self.tool, self.result_file_path[entry_name]).parse()
        self.ranked_mappings = get_ranked_mappings(self.mappings, self.unp_dir)
        end_cif = timer()
        logger.info(f'Total (from mmcif parsing to result parsing): {end_cif - start_cif} seconds.')

def run():
    parser = argparse.ArgumentParser(
        description="SIFTS mapping between a structure and a sequence."
    )

    parser.add_argument(
        "-i",
        "--cif-file",
        required=True,
        help="Base location for the PDBx/mmCIF file.",
    )
    parser.add_argument(
        "-od",
        "--output-dir",
        required=True,
        help="Base location for output SIFTS files.",
    )
    parser.add_argument(
        "-db",
        "--db-file",
        required=True,
        help="Location of the fasta(.gz) sequence database file.",
    )
    parser.add_argument(
        "-t",
        "--tool",
        default='mmseqs',
        help="Tool to use for the global alignment: blastp or mmseqs.",
    )
    parser.add_argument(
        "-ogm",
        "--out-global-mappings",
        required=False,
        help="Location of the csv file to save global mappings.",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use.",
    )    
    args = parser.parse_args()

    logger.info(vars(args))
    sifts_global_mappings = SiftsGlobalMappings(
        args.cif_file,
        args.output_dir,
        args.db_file,
        args.tool,
        threads=args.threads,
    )
    sifts_global_mappings.process()

if __name__ == "__main__":
    run()
