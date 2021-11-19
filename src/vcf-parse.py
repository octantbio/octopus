import re
import sys
import argparse
import itertools
from collections import Counter
from pathlib import Path

# catch broken pipe errors to allow ex) python foo.py ... | head
# see: https://stackoverflow.com/a/30091579
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

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

if __name__ == '__main__':
    description = "extract all barcodes vcf file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("fasta",
                        type=argparse.FileType('r'),
                        help="input fasta to parse for barcodes")
    parser.add_argument("vcf",
                        nargs="?",
                        type=argparse.FileType('r'),
                        help="vcf to parse (stdin if none)")
    args = parser.parse_args()

    # find all barcodes in the input fasta and drop their lens into a dict
    bcs = re.compile(r'N+')
    input_fasta = ((x[0], x[1].upper()) for x in fasta_reader(args.fasta))
    bc_dict = {
        ref: [len(x[0]) for x in re.finditer(bcs, seq)]
        for ref, seq in input_fasta
    }

    barcodes = []
    n_vars = 0
    for line in args.vcf:
        # use the reference to lookup barcode locations
        if line.startswith("##reference="):
            ref = Path(line).stem
            bc_lens = bc_dict[ref]
            expected_bcs = len(bc_lens)
        # ensure that len(variant) is in our lookup and the ref is all N's
        # this avoids a freebayes error where a mutation right next to the barcode would get looped in
        # we avoid checking positions as they can easily be off
        # pop the length from the lookup to avoid any false positives on the nth check
        # keep relevant portion of VCF file: [ref, pos, -, ref_seq, read_seq]
        if not line.startswith("#"):
            n_vars += 1
            vcf = line.rstrip().split('\t')[:5]
            if len(vcf[4]) in bc_lens and len(Counter(vcf[3])) == 1:
                barcodes.append((vcf[0], vcf[1], vcf[4]))
                bc_lens.remove(len(vcf[4]))

    # graceful exit for no barcode
    n_barcodes = len(barcodes)
    if n_barcodes == 0:
        print(f'{ref}\t1\tNA\tNA\t{n_vars}\t0\t{expected_bcs}',
              file=sys.stdout)
    # taking advantage of new f-strings
    for idx, tmp in enumerate(barcodes):
        ref, pos, bc = tmp
        print(
            f'{ref}\t{idx+1}\t{bc}\t{pos}\t{n_vars}\t{n_barcodes}\t{expected_bcs}',
            file=sys.stdout)
