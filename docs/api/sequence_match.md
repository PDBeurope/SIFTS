# Global mappings internals

Internal classes for building target databases and running sequence alignments
(BLAST or MMseqs2).  These are orchestrated by
[`SiftsSequenceMatch`](sifts_sequence_match.md).

## Abstract base classes

::: pdbe_sifts.sequence_match.database.ToolDatabase

::: pdbe_sifts.sequence_match.base_alignment_search.AlignmentSearch

## BLAST backend

::: pdbe_sifts.sequence_match.makeblastdb.MakeBlastDb

::: pdbe_sifts.sequence_match.blastp.BlastP

## MMseqs2 backend

::: pdbe_sifts.sequence_match.mmseqs_search.MmSearch

## Mappings parser

::: pdbe_sifts.sequence_match.sequence_match_parser.SequenceMatchParser

## Scoring helpers

::: pdbe_sifts.sequence_match.scoring_function_helper.get_tax_weight
