# ExportSIFTSTommCIF

Injects SIFTS mappings into mmCIF files by populating the `_pdbx_sifts_xref_db_segments` and `_pdbx_sifts_xref_db` categories, producing files compatible with official PDBe mmCIF releases.

::: pdbe_sifts.sifts_to_mmcif.main.ExportSIFTSTommCIF

## Supporting modules

### Mapping changes

::: pdbe_sifts.sifts_to_mmcif.delta_mappings.FindMappingChanges

### CSV / DB readers

::: pdbe_sifts.sifts_to_mmcif.read_sifts_csv.get_unp_segments

::: pdbe_sifts.sifts_to_mmcif.read_sifts_csv.get_unpres_mapping

### mmCIF table builder

::: pdbe_sifts.sifts_to_mmcif.comm_utils.get_xref_db

::: pdbe_sifts.sifts_to_mmcif.comm_utils.transpose_list

::: pdbe_sifts.sifts_to_mmcif.comm_utils.make_opt
