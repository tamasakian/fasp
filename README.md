# fasp

FASP: A FASTA Processor

## Install

```
pip3 install git+https://github.com/tamasakian/fasp.git
```

## Usage

```
python3 -m fasp <function> <args>
```

For example, you want to extract upstream region;

```
python3 -m fasp slice_sequence_by_upstream_region \
    input_filename \
    output_filename \
    output_id \
    output_description \
    strand \
    start \
    end \
    bp
```