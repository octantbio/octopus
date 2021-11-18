import re
import os
import sys
import argparse
import itertools

# catch broken pipe errors to allow ex) python foo.py ... | head
# see: https://stackoverflow.com/a/30091579
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

#===============================================================================

# silly over optimization to make a fast reverse compliment
# see: https://bioinformatics.stackexchange.com/q/3583
COMP = str.maketrans("ACTGacgt", "TGACtgca")


def rev_comp(seq):
    return seq.translate(COMP)[::-1]


#-------------------------------------------------------------------------------


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
    fa_iter = (x[1]
               for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read


#-------------------------------------------------------------------------------


def find_seq(pattern, seq):
    """
    Find any pattern on a plasmid (even taking into account sequences that span the junction)

    Input:
    ------
    pattern :: str
        pattern to search
    seq :: str
        sequence to search

    Output:
    -------
    loc :: int
        location of match

    Depends:
    --------
    itertools, re
    """
    search = [
        m.start(0) for m in re.finditer(pattern, seq + seq)
        if m.start() < len(seq)
    ]
    if not search:
        return -1
    else:
        return search[0]


#-------------------------------------------------------------------------------


def orient_seqs(pattern, fasta):
    """
    Re-orient sequences based on a pattern match

    Input:
    ------
    pattern :: str
        pattern to search
    fasta :: iter (name, seq)
        iterable that produces sequences

    Output:
    -------
    (good_seqs, bad_seqs)
    """
    good_seqs = []
    bad_seqs = []
    for name, seq in fasta:
        seq = seq.upper()
        ori_pos = find_seq(pattern, seq)
        if ori_pos != -1:
            good_seqs.append((name, seq[ori_pos:] + seq[:ori_pos]))
        else:
            # try again with the rc
            seq_rc = rev_comp(seq)
            ori_pos_2 = find_seq(pattern, seq_rc)
            if ori_pos_2 == -1:
                bad_seqs.append((name, seq))
            else:
                good_seqs.append(
                    (name, seq_rc[ori_pos_2:] + seq_rc[:ori_pos_2]))

    return (good_seqs, bad_seqs)


if __name__ == '__main__':
    description = "reorient all sequences in a fasta to start at a particular pattern"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("fasta",
                        nargs="+",
                        type=argparse.FileType('r'),
                        help="fasta's to process")
    parser.add_argument(
        '-r',
        '--records',
        type=int,
        dest='num_records',
        metavar='N',
        default='0',
        help='number of records in the fasta to process (default = all)')
    parser.add_argument('--ori',
                        type=str,
                        dest='ori_seq',
                        metavar='N',
                        default='GAAGATCCTTTGATT',
                        help='ori sequence (default = GAAGATCCTTTGATT)')
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        dest='verbose',
                        help='verbose logging')
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
        basename = os.path.basename(fasta.name).split('.')[0]
        records = itertools.islice(fasta_reader(fasta), num_records)
        good, bad = orient_seqs(args.ori_seq.upper(), records)

        # create the fasta header and include the file name
        for name, seq in good:
            print('>{}\n{}'.format(basename + "_" + name, seq),
                  file=sys.stdout)

        # send the bad reads to std.err if verbose
        if bad and args.verbose:
            for name, seq in bad:
                print('{}\t{}'.format(fasta.name, name), file=sys.stderr)
