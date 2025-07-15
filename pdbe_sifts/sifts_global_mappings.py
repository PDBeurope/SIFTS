#!/usr/bin/env python3

import argparse
import io
import csv
import json
from timeit import default_timer as timer
from pathlib import Path
import pickle

from pymmseqs.config.convertalis_config import ConvertAlisConfig

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif.entry import Entry
from pdbe_sifts.global_mappings.query_db import QueryDb
from pdbe_sifts.global_mappings.alignment_db import AlignDb
from pdbe_sifts.global_mappings.mmseqs_search import MmSearch
from pdbe_sifts.global_mappings.blastp import BlastP
from pdbe_sifts.global_mappings.alignment_result_parser import GlobMappingsParser
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
        self.fasta_files_path.mkdir(parents=True, exist_ok=True)
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
        self.queries = {}
        self.searches = {}
        self.aligns = {}
        self.convertalis = {}
        self.result_file_path = {}
        self.mappings = {}

    def create_mmseqs_query(self, id, fake_fasta):
        query_db_path = make_path(self.out_dir, id, 'queryDB', 'query.db')
        query_db = QueryDb(fake_fasta, query_db_path)
        query_db.run()
        self.queries[id] = query_db

    def mmseqs_search(self, id, query, target):
        now = get_date()
        result_db = make_path(self.out_dir, id, 'resultDB', 'result.db', now)
        tmp_fold = Path(str(self.out_dir / f'resultDB_{id}_{now}'/ 'tmp'))
        tmp_fold.parent.mkdir(parents=True, exist_ok=True)
        mm_search = MmSearch(query, target, result_db, tmp_fold)
        mm_search.run()
        self.searches[id] = mm_search
    
    def blastp_search(self, id, fasta_path):
        output_path = make_path(self.out_dir, id, 'blastp', 'blast_result.json')
        blastp_search = BlastP(fasta_path, self.db_file, output_path)
        blastp_search.run()
        self.searches[id] = blastp_search
        self.result_file_path[id] = output_path

    def mmseqs_align(self, id, query, target, result):
        align_db_path = make_path(self.out_dir, id, 'alignDB', 'align.db')
        align_db = AlignDb(result, align_db_path, query, target, v=0)
        align_db.run()
        self.aligns[id] = align_db
    
    def mmseqs_convertalis(self, id, query, target, align):
        convertalis_path = make_path(self.out_dir, id, 'hits', f'hits_{id}.tsv')
        format_string = 'query,target,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,qaln,taln,qlen,taxid'
        convertalis_pro = ConvertAlisConfig(query, target, align, convertalis_path, format_output=format_string)
        convertalis_pro.run()
        self.convertalis[id] = convertalis_pro
        self.result_file_path[id] = convertalis_path

    def search(self, id, tmp_fasta_path):
        match self.tool:
            case 'mmseqs':
                self.create_mmseqs_query(id, tmp_fasta_path)
                self.mmseqs_search(id, self.queries[id].query_path, self.db_file)
                self.mmseqs_align(id, self.queries[id].query_path, self.db_file, self.searches[id].output_path)
                self.mmseqs_convertalis(id, self.queries[id].query_path, self.db_file, self.aligns[id].output_path)
            case 'blastp':
                self.blastp_search(id, tmp_fasta_path)
            case _:
                raise ValueError(f"Unsupported tool: {self.tool}")

    def process(self):
        logger.info("Processing [%s]" % self.cif_file)
        start = timer()
        entry_name = Path(self.cif_file).stem
        self.entry = Entry(entry_name, self.cif_file)
        entity_seq_tax = self.entry.get_entity_seq_tax()
        for ent, seq_tax in entity_seq_tax.items():
            fake_header = f'>pdb|{entry_name}-{ent}|None'
            fake_seq = seq_tax[0]
            tax_id = seq_tax[1]
            fake_content = f'{fake_header}\n{fake_seq}'
            tmp_fasta_path = self.fasta_files_path / f'tmp_{entry_name}_{ent}.fasta'
            with open(tmp_fasta_path, 'w') as f:
                f.write(fake_content)
            id = f'{entry_name}-{ent}'
            self.search(id, tmp_fasta_path)
            self.mappings[id] = GlobMappingsParser(self.tool, id, self.result_file_path[id], tax_id).parse()
            ranked_mappings = get_ranked_mappings(self.mappings[id])
            print(ranked_mappings)
            end = timer()
            logger.info(f'Total (from mmcif parsing to result parsing): {end-start} seconds.')

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

    args = parser.parse_args()

    logger.info(vars(args))
    sifts_global_mappings = SiftsGlobalMappings(
        args.cif_file,
        args.output_dir,
        args.db_file,
        args.tool,
    )
    sifts_global_mappings.process()

if __name__ == "__main__":
    run()
