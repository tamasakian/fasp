#!/usr/bin/env python3

import argparse
import sys
from fasp import fasta
from fasp import fastp
from fasp import fastn

def main():
    args = parse_args()
    function = functions[args.function]
    try:
        function(*args.args)
    except TypeError as e:
        print(e)
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('function', choices=functions)
    parser.add_argument('args', type=str, nargs='*')
    args = parser.parse_args()
    return args

functions = {
    "rename_header": fasta.rename_header, 
    "prefix_to_sequence_ids": fasta.prefix_to_sequence_ids, 
    "slice_records_by_exact_ids": fasta.slice_records_by_exact_ids, 
    "slice_records_by_partial_ids": fasta.slice_records_by_partial_ids, 
    "rename_headers_feature": fastn.rename_headers_feature
}

if __name__ == "__main__":
    main()
