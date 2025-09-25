import csv
import json

IDENTITY_CUTOFF = 0.9

class GlobMappingsParser:
    def __init__(self, format, result_file_path):
        self.format = format
        self.result_file_path = result_file_path
        self.mappings = None

    def parse(self):
        if self.format == 'mmseqs':
            self._parse_mmseqs()
        elif self.format == 'blastp':
            self._parse_blastp()
        else:
            raise ValueError(f"Unsupported format: {format}")
        return self.mappings

    def _build_mapping_dict(
        self,
        entry,
        entity,
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
        query_tax_id,
        target_tax_id,
    ):
        return {
            'entry': entry,
            'entity': entity,
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
            'query_tax_id': query_tax_id,
        }

    def _parse_mmseqs(self):
        data = {}
        with open(self.result_file_path) as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            next(reader)
            for row in reader:
                # >pdb|entry-entity|OX=taxid
                # 'query,target,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,qaln,taln,qlen,taxid,qheader,fident,qcov'
                header_list = row[14].split('|')
                entry, entity = header_list[1].split('-')
                taxid = header_list[-1].split('=')[-1]
                if entry not in data:
                    data[entry]= {}
                if entity not in data[entry]:
                    data[entry][entity] = []
                taxid = int(taxid) if taxid not in (None, "None", "") else -1
                identity_crt = float(row[15])
                if identity_crt >= IDENTITY_CUTOFF:
                    data[entry][entity].append(self._build_mapping_dict(
                        entry=entry,
                        entity=entity,
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
                        coverage=float(row[16]),
                        query_aligned=row[10],
                        target_aligned=row[11],
                        target_tax_id= int(row[13]),
                        query_tax_id=taxid
                    ))
        self.mappings = data

    def _parse_blastp(self):
        data = {}
        with open(self.result_file_path) as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                # >pdb|entry-entity|OX=taxid
                #  qseqid sseqid length mismatch qstart qend sstart send evalue bitscore qseq sseq qlen staxid pident qcovs
                header_list = row[0].split('|')
                entry, entity = header_list[1].split('-')
                if entry not in data:
                    data[entry]= {}
                if entity not in data[entry]:
                    data[entry][entity] = []
                taxid = int(header_list[-1].split('=')[-1])
                identity_crt = float(row[14])/100
                if identity_crt >= IDENTITY_CUTOFF:
                    data[entry][entity].append(self._build_mapping_dict(
                        entry=entry,
                        entity=entity,
                        accession=row[1].split('|')[1],# unp|acc|name
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
                        coverage=float(row[15])/100,
                        query_aligned=row[10],
                        target_aligned=row[11],
                        target_tax_id= int(row[13]),
                        query_tax_id=taxid
                    ))
        self.mappings = data