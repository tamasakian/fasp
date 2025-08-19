# fasp

FASP: A FASTA Processor

## Install

```
pip3 install git+https://github.com/tamasakian/fasp.git
```

If you want to install in editable mode;
```
pip3 install -e git+https://github.com/tamasakian/fasp.git#egg=fasp
```

## Syntax

```
python3 -m fasp <function> <args>
```

## Functions
`slice_tail`: Slice the tail of a FASTA record.

```
python3 -m fasp slice_tail <input_fasta> <output_fasta> <length>
```
- `input_fasta`: Input FASTA file.
- `output_fasta`: Output FASTA file.
- `length`: Length of the tail to slice.