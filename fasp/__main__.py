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
    "slice_records_by_exact_ids": fasta.slice_records_by_exact_ids, 
    "slice_records_by_partial_ids": fasta.slice_records_by_partial_ids
}

if __name__ == "__main__":
    main()
