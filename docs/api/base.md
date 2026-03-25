# Base utilities

Low-level helpers shared across the pipeline: path resolution, CPU detection, and all custom exceptions.

## Path helpers

::: pdbe_sifts.base.paths.uniprot_cache_dir

::: pdbe_sifts.base.paths.ccd_cache_path

## CPU helpers

::: pdbe_sifts.base.utils.get_allocated_cpus

::: pdbe_sifts.base.utils.get_cpu_count

## Exceptions

::: pdbe_sifts.base.exceptions.ObsoleteUniProtError

::: pdbe_sifts.base.exceptions.AccessionNotFound

::: pdbe_sifts.base.exceptions.NotAPolyPeptide

::: pdbe_sifts.base.exceptions.NoSegmentsError

::: pdbe_sifts.base.exceptions.SplitAccessionError

::: pdbe_sifts.base.exceptions.BatchRunException

::: pdbe_sifts.base.exceptions.EntryFailedException

::: pdbe_sifts.base.exceptions.ReleaseCheckFailedException

::: pdbe_sifts.base.exceptions.ProcessFailedError

::: pdbe_sifts.base.exceptions.EntryTimedOutException

## Logging

::: pdbe_sifts.base.log.SensitiveFormatter

## Argparse action

::: pdbe_sifts.base.utils.SiftsAction
