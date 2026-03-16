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
