# mmCIF parsing

Low-level classes for reading and representing PDB mmCIF files.
`Entry` is the main entry point; it parses a CIF file and exposes `Chain` objects,
each backed by an `Entity` that holds residue-level data.

## Entry

::: pdbe_sifts.mmcif.entry.Entry

## Chain

::: pdbe_sifts.mmcif.chain.Chain

## Entity

::: pdbe_sifts.mmcif.entity.Entity

## Residue

::: pdbe_sifts.mmcif.residue.Residue

## ChemCompMapping

::: pdbe_sifts.mmcif.chem_comp.ChemCompMapping

## mmCIF (low-level parser)

::: pdbe_sifts.mmcif.mmcif_helper.mmCIF
