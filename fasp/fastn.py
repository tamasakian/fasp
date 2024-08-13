#!/usr/bin/env python3

"""Library for processing nucleotide FASTA files.

Functions
---------
rename_headers_feature: 
    Rename headers feature values.

"""

import re
from Bio import SeqIO

def rename_headers_feature(input_filename: str, output_filename: str, feature: str) -> None:
    """Rename headers feature values.

    Args:
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    feature : str
        Feature to rename headers. 
        e.g., locus_tag, protein_id.

    """
    pattern = re.compile(r"\[" + re.escape(feature) + r"=([^\]]+)]")

    def extract_feature_value(text: str, default: str = "") -> str:
        """Extract feature value from text if it matches the pattern.

        Args
        ----
        text : str
            Text to extract feature value.
        default : str
            Default value which is returned if the feature is not found. 
        
        Returns
        -------
        match.group(1) : str
            Feature value if the feature is found.
            If this is not found, default is returned instead.

        """
        match = pattern.search(text)
        return match.group(1) if match else default

    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            feature_value = (extract_feature_value(record.description) or
                             extract_feature_value(record.name) or
                             extract_feature_value(record.id, default=record.id))
            record.id = feature_value
            record.name = ""
            record.description = ""  
            SeqIO.write(record, output_handle, "fasta")