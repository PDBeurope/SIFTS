# Global mappings internals

Internal classes for building target databases and running sequence alignments
(BLAST or MMseqs2).  These are orchestrated by
[`SiftsGlobalMappings`](sifts_global_mappings.md).

## Abstract base classes

::: pdbe_sifts.global_mappings.database.ToolDatabase

::: pdbe_sifts.global_mappings.base_alignment_search.AlignmentSearch

## BLAST backend

::: pdbe_sifts.global_mappings.makeblastdb.MakeBlastDb

::: pdbe_sifts.global_mappings.blastp.BlastP

## MMseqs2 backend

::: pdbe_sifts.global_mappings.mmseqs_search.MmSearch

## Mappings parser

::: pdbe_sifts.global_mappings.global_mappings_parser.GlobMappingsParser

## Scoring helpers

::: pdbe_sifts.global_mappings.scoring_function_helper.get_tax_weight
