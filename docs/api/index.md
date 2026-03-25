# API Reference

All pipeline stages can be driven directly from Python without going through the CLI.

## Public classes

| Class | Module | Purpose |
|-------|--------|---------|
| [`SiftsGlobalMappings`](sifts_global_mappings.md) | `pdbe_sifts.sifts_global_mappings` | Alignment and scoring pipeline |
| [`SiftsAlign`](sifts_segments_generation.md) | `pdbe_sifts.sifts_segments_generation` | Per-entry segment and residue mapping |
| [`FastaBuilder`](sifts_fasta_builder.md) | `pdbe_sifts.sifts_fasta_builder` | Extract sequences from mmCIF files into a FASTA |
| [`TargetDb`](target_database.md) | `pdbe_sifts.global_mappings.target_database` | Build a MMseqs2 or BLAST reference database |
| [`SiftsDB`](sifts_db_wrapper.md) | `pdbe_sifts.database.sifts_db_wrapper` | DuckDB schema management and bulk loader |
| [`ExportSIFTSTommCIF`](sifts_to_mmcif.md) | `pdbe_sifts.sifts_to_mmcif.main` | Write SIFTS-annotated mmCIF files |
| [`UNP`](unp.md) | `pdbe_sifts.unp.unp` | UniProt REST client and isoform resolver |
| [`Entry`](mmcif.md) | `pdbe_sifts.mmcif.entry` | Parse a PDB mmCIF file into chain objects |
| [`Chain`](mmcif.md) | `pdbe_sifts.mmcif.chain` | Per-chain residue maps and segment data |
| [`Entity`](mmcif.md) | `pdbe_sifts.mmcif.entity` | Per-entity sequence and isoform alignment |
| [`AlignmentSearch`](global_mappings.md) | `pdbe_sifts.global_mappings.base_alignment_search` | Abstract base for BLAST/MMseqs2 backends |
| [`BlastP`](global_mappings.md) | `pdbe_sifts.global_mappings.blastp` | BLAST alignment backend |
| [`MmSearch`](global_mappings.md) | `pdbe_sifts.global_mappings.mmseqs_search` | MMseqs2 alignment backend |
| [`GlobMappingsParser`](global_mappings.md) | `pdbe_sifts.global_mappings.global_mappings_parser` | Parse and score raw alignment hits |
| [`EntryMapping`](segments_generation.md) | `pdbe_sifts.segments_generation.alignment.helper` | Per-entry segment/residue alignment orchestrator |
| [`ConnectivityCheck`](segments_generation.md) | `pdbe_sifts.segments_generation.connectivity.process_connectivity` | Check and refine segment connectivity |
| [`BatchSegmentGeneration`](batch.md) | `pdbe_sifts.batch_segment_generation` | Run segment generation over many entries |
| [`TaxonomyFix`](batch.md) | `pdbe_sifts.taxonomy_fix_pkl` | Resolve taxonomy IDs from bundled pickle |

## Typical usage pattern

```python
from pdbe_sifts.global_mappings.target_database import TargetDb
from pdbe_sifts.sifts_global_mappings import SiftsGlobalMappings
from pdbe_sifts.sifts_segments_generation import SiftsAlign
import duckdb
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB

# 1. Build reference DB (once)
TargetDb("uniprot_sprot.fasta", "./db", "taxonomy.tsv").run()

# 2. Global mappings
SiftsGlobalMappings("1abc.cif", "./results", "./db/target_db").process()

# 3. Segment generation
sa = SiftsAlign("1abc.cif", "./segments", db_conn_str="./results/hits.duckdb")
sa.process_entry("1abc")

# 4. Load into DuckDB
conn = duckdb.connect("./results/hits.duckdb")
SiftsDB(conn).bulk_load_from_entries("./segments")
conn.close()
```
