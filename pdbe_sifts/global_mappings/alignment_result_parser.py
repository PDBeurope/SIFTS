import csv
import json

from pdbe_sifts.base.utils import get_mismatches, get_identity, get_coverage


IDENTITY_CUTOFF = 0.9

class GlobMappingsParser:
    def __init__(self, format, id, entry, entity, result_file_path, query_tax_id):
        self.id = id
        self.entry = entry
        self.entity = entity
        self.format = format
        self.result_file_path = result_file_path
        self.mappings = None
        self.query_tax_id = query_tax_id

    def parse(self):
        if self.format == 'mmseqs':
            self._parse_mmseqs(self.id)
        elif self.format == 'blastp':
            self._parse_blastp(self.id)
        else:
            raise ValueError(f"Unsupported format: {format}")
        return self.mappings

    def _build_mapping_dict(
        self,
        entry_entity,
        accession,
        alignment_len,
        query_len,
        mismatch,
        identity,
        coverage,
        query_start,
        query_end,
        target_start,
        target_end,
        evalue,
        bit_score,
        query_aligned,
        target_aligned,
        # query_tax_id,
        target_tax_id = None,
    ):
        return {
            'entry-entity': entry_entity,
            'entry': self.entry,
            'entity': self.entity,
            'accession': accession,
            'alignment_len': alignment_len,
            'query_len': query_len,
            'mismatch': mismatch,
            'identity': identity,
            'coverage': coverage,
            'query_start': query_start,
            'query_end': query_end,
            'target_start': target_start,
            'target_end': target_end,
            'evalue': evalue,
            'bit_score': bit_score,
            'query_aligned': query_aligned,
            'target_aligned': target_aligned,
            'target_tax_id': target_tax_id,
            'query_tax_id': self.query_tax_id,
        }

    def _parse_mmseqs(self, id):
        data = []
        with open(self.result_file_path) as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                identity_crt = get_identity(row[10], row[11])
                if identity_crt >= IDENTITY_CUTOFF:
                    data.append(self._build_mapping_dict(
                        entry_entity=row[0],
                        accession=row[1],
                        alignment_len=int(row[2]),
                        query_len=int(row[12]),
                        mismatch=int(row[3]),
                        query_start=int(row[4]),
                        query_end=int(row[5]),
                        target_start=int(row[6]),
                        target_end=int(row[7]),
                        evalue=float(row[8]),
                        bit_score=float(row[9]),
                        identity=identity_crt,
                        coverage=get_coverage(int(row[4]), int(row[5]), int(row[12])),
                        query_aligned=row[10],
                        target_aligned=row[11],
                        target_tax_id= int(row[13]),
                    ))
        self.mappings = data

    def _parse_blastp(self, id):
        data = []
        with open(self.result_file_path) as jfile:
            jfile_dict = json.load(jfile)

        results = jfile_dict['BlastOutput2'][0]['report']['results']['search']
        qlen = results['query_len']
        entry_entity = results['query_title'].split('|')[1]

        for hit in results['hits']:
            hsp = hit['hsps'][0]
            identity_crt = get_identity(hsp['qseq'], hsp['hseq'])
            if identity_crt >= IDENTITY_CUTOFF:
                accession = hit['description'][0]['accession']
                target_tax_id = hit['description'][0]['taxid']
                data.append(self._build_mapping_dict(
                    entry_entity=entry_entity,
                    accession=accession,
                    alignment_len=hsp['align_len'],
                    query_len=qlen,
                    mismatch=get_mismatches(hsp['qseq'], hsp['hseq']),
                    identity=identity_crt,
                    coverage=get_coverage(hsp['query_from'], hsp['query_to'], qlen),
                    query_start=hsp['query_from'],
                    query_end=hsp['query_to'],
                    target_start=hsp['hit_from'],
                    target_end=hsp['hit_to'],
                    evalue=hsp['evalue'],
                    bit_score=hsp['bit_score'],
                    query_aligned=hsp['qseq'],
                    target_aligned=hsp['hseq'],
                    target_tax_id=target_tax_id,
                ))
        self.mappings = data