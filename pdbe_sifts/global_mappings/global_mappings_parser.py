import csv
import multiprocessing as mp
import pickle
from pathlib import Path
from pdbe_sifts.global_mappings.scoring_function import get_ranked_mappings

IDENTITY_CUTOFF = 0.9

def write_mapping(loc, data):
    """Serialize parsed or ranked mapping to a pickle file."""
    with open(str(loc), 'wb') as f:
        pickle.dump(data, f)


def _parse_blastp_rows(rows, entry):
    """Parse BLASTP tabular rows (16 columns) into the unified data structure."""
    data = {}
    for row in rows:
        qseqid = row[0]
        header_list = qseqid.split('|')
        _, entity = header_list[1].split('-')

        taxid = header_list[-1].split('=')[-1]
        taxid = int(taxid) if taxid.isdigit() else -1

        identity_crt = float(row[14]) / 100.0  # pident (% → fraction)
        coverage = float(row[15]) / 100.0      # qcovs (% → fraction)
        if entry not in data:
            data[entry] = {}
        if entity not in data[entry]:
            data[entry][entity] = []
        if identity_crt >= IDENTITY_CUTOFF:
            data[entry][entity].append({
                'entry': entry,
                'entity': entity,
                'accession': row[1].split('|')[1],# unp|acc|name
                'alignment_len': int(row[2]),
                'query_len': int(row[12]),
                'mismatch': int(row[3]),
                'query_start': int(row[4]),
                'query_end': int(row[5]),
                'target_start': int(row[6]),
                'target_end': int(row[7]),
                'evalue': float(row[8]),
                'bit_score': float(row[9]),
                'identity': identity_crt,
                'coverage': coverage,
                'query_aligned': row[10],
                'target_aligned': row[11],
                'target_tax_id': int(row[13]),
                'query_tax_id': taxid,
            })
    return data

def _parse_mmseqs_rows(rows, entry):
    """Parse MMseqs2 tabular rows (17 columns) into the unified data structure."""
    data = {}
    for row in rows:
        header_list = row[14].split('|')
        _, entity = header_list[1].split('-')
        taxid = header_list[-1].split('=')[-1]
        taxid = int(taxid) if taxid.isdigit() else -1
        identity_crt = float(row[15])
        coverage = float(row[16])
        if entry not in data:
            data[entry] = {}
        if entity not in data[entry]:
            data[entry][entity] = []
        if identity_crt >= IDENTITY_CUTOFF:
            data[entry][entity].append({
                'entry': entry,
                'entity': entity,
                'accession': row[1],
                'alignment_len': int(row[2]),
                'query_len': int(row[12]),
                'mismatch': int(row[3]),
                'query_start': int(row[4]),
                'query_end': int(row[5]),
                'target_start': int(row[6]),
                'target_end': int(row[7]),
                'evalue': float(row[8]),
                'bit_score': float(row[9]),
                'identity': identity_crt,
                'coverage': coverage,
                'query_aligned': row[10],
                'target_aligned': row[11],
                'target_tax_id': int(row[13]),
                'query_tax_id': taxid,
            })
    return data

def _process_entry(entry, entity, rows, unp_dir, out_dir_parsed, out_dir_ranked):
    """Worker function to process all rows for one (entry, entity) pair."""
    if not rows:
        return

    n_cols = len(rows[0])
    if n_cols == 16:
        data = _parse_blastp_rows(rows, entry)
    elif n_cols == 17:
        data = _parse_mmseqs_rows(rows, entry)
    else:
        return  # unsupported format

    # Write parsed pickle
    parsed_path = out_dir_parsed / f"{entry}_{entity}.pkl"
    write_mapping(parsed_path, data)

    # Write ranked pickle
    ranked = get_ranked_mappings(data, unp_dir)
    if ranked:
        ranked_path = out_dir_ranked / f"{entry}_{entity}.pkl"
        write_mapping(ranked_path, ranked)


def _process_entry_wrapper(args):
    """Unpack tuple args for pool.imap_unordered."""
    return _process_entry(*args)


def _entry_generator(result_file_path):
    """
    Stream the TSV file and yield (entry, entity, [rows]) for each unique (entry, entity) pair.

    The output files from MMseqs2/BLASTP may not be globally ordered by entry
    when run with multiple threads, but all rows for a given (entry, entity)
    are guaranteed to be contiguous.
    """
    with open(result_file_path) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        current_entry, current_entity = None, None
        buffer = []

        for row in reader:
            # header is in column 0 for BLASTP, column 14 for MMseqs2
            header = row[14] if len(row) == 17 else row[0]
            entry, entity = header.split('|')[1].split('-')

            # Initialize on first line
            if current_entry is None and current_entity is None:
                current_entry, current_entity = entry, entity

            # When we encounter a new (entry, entity) pair → yield the previous buffer
            if (entry != current_entry) or (entity != current_entity):
                yield current_entry, current_entity, buffer
                buffer = []
                current_entry, current_entity = entry, entity

            buffer.append(row)

        # yield the last group
        if buffer:
            yield current_entry, current_entity, buffer



class GlobMappingsParser:
    def __init__(self, format, result_file_path, out_dir, unp_dir, max_workers=None):
        self.format = format
        self.result_file_path = result_file_path
        self.mappings = None
        self.out_dir = Path(out_dir)
        self.unp_dir = unp_dir
        self.out_dir_parsed = self.out_dir / 'parsed_hits'
        self.out_dir_ranked = self.out_dir / 'ranked_hits'
        self.out_dir_parsed.mkdir(parents=True, exist_ok=True)
        self.out_dir_ranked.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers or mp.cpu_count()

    def parse(self):
        """Dispatch to the correct parsing method based on file format."""
        if self.format in ('mmseqs', 'blastp'):
            self._parse_generic()
        else:
            raise ValueError(f"Unsupported format: {self.format}")
        return self.mappings

    def _parse_generic(self):
        """Generic parser for MMseqs2 and BLASTP outputs (grouped by entry + entity)."""
        tasks = (
            (entry, entity, rows, self.unp_dir, self.out_dir_parsed, self.out_dir_ranked)
            for entry, entity, rows in _entry_generator(self.result_file_path)
        )

        with mp.Pool(processes=self.max_workers) as pool:
            for _ in pool.imap_unordered(_process_entry_wrapper, tasks):
                pass
