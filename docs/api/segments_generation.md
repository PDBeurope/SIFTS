# Segments generation

Modules that convert global alignment hits into per-residue SIFTS mappings
and write the segment / residue CSV files consumed by the rest of the pipeline.

## Alignment helpers

### Data types

::: pdbe_sifts.segments_generation.alignment.helper.SMapping

::: pdbe_sifts.segments_generation.alignment.helper.SiftsMode

### CustomSequenceAccession

::: pdbe_sifts.segments_generation.alignment.helper.CustomSequenceAccession

### EntryMapping

::: pdbe_sifts.segments_generation.alignment.helper.EntryMapping

### Alignment functions

::: pdbe_sifts.segments_generation.alignment.do_alignment_lalign36

::: pdbe_sifts.segments_generation.alignment.annotate_alignment

::: pdbe_sifts.segments_generation.alignment.remove_range_alignment

::: pdbe_sifts.segments_generation.alignment.get_conflicts

::: pdbe_sifts.segments_generation.alignment.get_identity

::: pdbe_sifts.segments_generation.alignment.get_align_chunk

### Utility functions

::: pdbe_sifts.segments_generation.alignment.helper.get_accession

::: pdbe_sifts.segments_generation.alignment.helper.fmt_ranges

::: pdbe_sifts.segments_generation.alignment.helper.overlapping

::: pdbe_sifts.segments_generation.alignment.helper.out_of_order

::: pdbe_sifts.segments_generation.alignment.helper.unroll_map

::: pdbe_sifts.segments_generation.alignment.helper.large_gap

## Connectivity

::: pdbe_sifts.segments_generation.connectivity.ccd_parser.CcdFile

::: pdbe_sifts.segments_generation.connectivity.process_connectivity.ConnectivityCheck

## CSV generation

### Data types

::: pdbe_sifts.segments_generation.generate_xref_csv.XRefSegment

::: pdbe_sifts.segments_generation.generate_xref_csv.XRefResidue

### Functions

::: pdbe_sifts.segments_generation.generate_xref_csv.insert_mappings

::: pdbe_sifts.segments_generation.generate_xref_csv.write_seg_csv

::: pdbe_sifts.segments_generation.generate_xref_csv.write_res_csv

::: pdbe_sifts.segments_generation.generate_xref_csv.insert_residues

## Mapping queries

::: pdbe_sifts.segments_generation.get_list_of_mappings.get_selection_queries

::: pdbe_sifts.segments_generation.get_list_of_mappings.get_curated_db_mappings
