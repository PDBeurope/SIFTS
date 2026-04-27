"""
Author: Preeti Choudhary
Desription: Parse the ccd cif file and generate three_to_one_letter_mapping.csv file
"""

import gemmi
import argparse
import csv
from pathlib import Path
from pdbe_sifts.base import utils



def parse_ccd_file(ccd_cif_file, ccd_mapping_file):
    
    """
    Description:Parse the PDBe CCD cif file and generate three to one letter mapping file
    Args: Input - ccd cif file (components.cif - can be downloaded from PDBe FTP)
    Output- three_to_one_letter_mapping.csv
    """
    # Open output three_to_one_letter_mapping.csv file for writing
    
    print(f"Writing the three_to_one_mapping.csv: {ccd_mapping_file}")
    with open(ccd_mapping_file, "w", newline="") as f:
        writer = csv.writer(f)

        # Read PDBe CCD CIF file
        # make sure file path is a string as gemmi does not accept a Path object;
        doc = gemmi.cif.read_file(str(ccd_cif_file))

        for block in doc:
            # get the chemp_comp category
            chem_comp_cat = block.get_mmcif_category("_chem_comp")

            if chem_comp_cat:
                # get the one/three letter code for a given ccd
                one_letter_list = chem_comp_cat.get("one_letter_code")
                three_letter_list = chem_comp_cat.get("three_letter_code")

                if one_letter_list and three_letter_list:
                    one = one_letter_list[0]
                    three = three_letter_list[0]

                    # If missing values i.e ("?" or ".") in use X as one letter code                  
                    if one is None :
                        one = "X"
                    writer.writerow([three,one])
        

# Path of the CCD cif file from PDBe FTP area
DEFAULT_CCD_URL='https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/ccd/components.cif'


def main():
    parser = argparse.ArgumentParser(
        description="Generate three_to_one_letter_mapping.csv from CCD CIF file"
    )

    parser.add_argument(
        "--ccd_file",
        help="Path to CCD CIF file (if not provided, it will be downloaded)",
    )

    parser.add_argument(
        "--output_path",
        required=True,
        help="Path to output three_to_one_mapping.csv file",
    )

    parser.add_argument(
        "--ccd_url",
        default=DEFAULT_CCD_URL,
        help="URL to download CCD file if --ccd_file is not provided",
    )

    args = parser.parse_args()
    

    # Decide where CCD file comes from
    if args.ccd_file:
        ccd_file = Path(args.ccd_file)
    else:
        # Default download location - current directory from where code is being run
        base_dir = Path.cwd()
        ccd_file = Path(base_dir,"components.cif")
        if not ccd_file.exists():
            utils.download_file_from_url(args.ccd_url, ccd_file)
        else:
            print(f"Using existing ccd file: {ccd_file}")
  
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = Path(output_path, "three_to_one_letter_mapping.csv")
    
    parse_ccd_file(ccd_file, output_file)


if __name__ == "__main__":
    main()