import re
import sys
import argparse
import itertools

from pathlib import Path

# catch broken pipe errors to allow ex) python foo.py ... | head
# see: https://stackoverflow.com/a/30091579
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence
    Input:
    ------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.
    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]
    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read

#===============================================================================

if __name__ == '__main__':
    description = "reorient all sequences in a fasta to start at a particular pattern"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("fasta",
                        nargs="+",
                        type=argparse.FileType('r'),
                        help="fasta's to process")
    parser.add_argument('-r',
                        '--records',
                        type=int,
                        dest='num_records',
                        metavar='N',
                        default='0',
                        help='number of records in the fasta to process (default = all)')
    parser.add_argument('--no-flat',
                        dest='no_flat',
                        action='store_true',
                        help='just output the file'
                        )
    args = parser.parse_args()

    if args.num_records < 0:
        raise ValueError('Number of records must be >0!!')
    elif args.num_records == 0:
        num_records = None
    else:
        num_records = args.num_records

    #---------------------------------------------------------------------------

    # re-orient the reference plasmids
    for fasta in args.fasta:
        basename = Path(fasta.name).stem
        records = itertools.islice(fasta_reader(fasta), num_records)

        # include the well in the record line
        # also concatenate sequences together to enable mapping across the plasmid junction
        for name, seq in records:
            seq = seq.upper()
            if args.no_flat:
                print(f'>{basename + "_" + name}\n{seq}', file=sys.stdout)
            else:
                print(f'>{basename + "_" + name}\n{seq+seq}', file=sys.stdout)

