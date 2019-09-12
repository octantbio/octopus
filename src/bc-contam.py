import re
import sys
import shlex
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool,cpu_count
from collections import Counter, defaultdict

# catch broken pipe errors to allow ex) python foo.py ... | head
# see: https://stackoverflow.com/a/30091579
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def proc_tview(pack):
    bam, ref, bc, loc = pack
    tview_cmd = f'samtools tview -d T -p {ref}:{loc}-{int(loc) + len(bc)} {bam}'
    tview = subprocess.run(
            shlex.split(tview_cmd),
            stdout=subprocess.PIPE
            )

    # normalize text, ignore tview header (first 3 lines),
    # trim barcode, only take full barcodes
    bcs_reads = tview.stdout.decode('utf-8').rstrip().upper().splitlines()[3:]
    bcs_reads = (x[:len(bc)] for x in bcs_reads)
    bcs_reads = Counter(x for x in bcs_reads if not re.search(r'[^ATGC]', x))
    coverage = sum(bcs_reads.values())

    # collapse bcs with starcode
    starcode_cmd = 'starcode --sphere --dist 1'
    starcode = subprocess.run(
            shlex.split(starcode_cmd),
            input='\n'.join(f'{bc}\t{reads}' for bc, reads in bcs_reads.items()).encode('utf-8'),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL
            )

    # process the output (split on new lines then tab)
    star_out = starcode.stdout.decode('utf-8').rstrip().splitlines()
    star_out = (x.split() for x in star_out)
    bcs_collapse = Counter({bc:reads for bc,reads in star_out})

    return((bam, (coverage, bcs_collapse)))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    description = "Report all reads at the highest depth barcode to check for contamination"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input",
                        nargs="?",
                        type=argparse.FileType('r'),
                        help="tidy list of bam file, what it aligned to, a barcode, and its location")
    parser.add_argument('-j',
                        '--proc',
                        dest='proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, cpu_count()+1),
                        help='number of processors (default=1, max={})'.format(cpu_count()))
    args = parser.parse_args()

    # build input
    tmp_input = (x.split()[:4] for x in args.input)
    to_run = [x for x in tmp_input if x[2] != 'NA']

    # 1) run samtools tview on each entry
    # 2) collapse all barcodes for a plasmid into well:[(coverage, Counter(bcs)), ...]
    well_dict = defaultdict(list)
    with Pool(args.proc) as pool:
        tview_out = pool.imap_unordered(proc_tview, to_run, chunksize=32)
        for bam, bc_collapse in tview_out:
            well_dict[bam].append(bc_collapse)

    # 3) find barcode with highest coverage
    # 4) output all barcodes with their reads at that site
    # recall bc_collapse = [(cov, Counter(reads)), (cov, Counter(reads)), ...]
    for bam, bc_collapse in well_dict.items():
        max_cov = max(bc_collapse, key=lambda x: x[0])
        for bc, reads in max_cov[1].items():
            print(f'{bam}\t{bc}\t{reads}', file=sys.stdout)
