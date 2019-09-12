import sys
import argparse
import itertools
import multiprocessing
from signal import signal, SIGPIPE, SIG_DFL
from collections import defaultdict

# external deps
import pygsheets
import Levenshtein

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

# authorize gsheets
gsheet = pygsheets.authorize(
        client_secret='/data/projects/platform/src/pygsheets.json',
        credentials_directory='/data/projects/platform/src'
        )

#-------------------------------------------------------------------------------

# silly over optimization to make a fast reverse compliment
# see: https://bioinformatics.stackexchange.com/q/3583
COMP = str.maketrans("ACTGacgt", "TGACtgca")
def rev_comp(seq):
    return seq.translate(COMP)[::-1]

#-------------------------------------------------------------------------------

# compare barcode to list of barcodes
def test_bc(pack):
    """
    Compare a single barcode to an input list and output any conflicts. Note we
    packed all of the relevant variables together to make mapping easier

    Input:
    ------
    bc :: str
        barcode to test
    dist :: int
        levenshtein distance cutoff
    in_list :: (db_name, barcode)
        a flattened list of barcodes to test with index and database name
    rest :: str
        relevant columns to identify what well the barcode came from

    Output:
    -------
    {rest:[(barcode, dist, db_name, barcode), ...]
        offending barcode, its distance, and what other barcodes it conflict with
        in a defaultdict
    """
    bc, dist, in_list, rest = pack
    conflict = defaultdict(list)
    for db_name, test_bc in in_list:
        this_dist = Levenshtein.distance(bc, test_bc)
        if this_dist <= dist:
            conflict[rest].append((bc, str(this_dist), db_name, test_bc))
    if not conflict:
        conflict[rest].append((bc, 'NA', 'NA', 'NA'))
    return(conflict)

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Check barcode(s) against our known databases. Operates on the last column of the input')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to file (or stdin if none)')
    parser.add_argument('-d',
                        '--dist',
                        dest='dist',
                        type=int,
                        default=3,
                        metavar='N',
                        help='min distance to report (default=3)')
    parser.add_argument('-f',
                        '--format',
                        choices=['tidy', 'quiet'],
                        default='quiet',
                        dest='format',
                        help='format output - quiet reports TRUE/FALSE')
    parser.add_argument('-t', '-j',
                        '--threads',
                        dest='proc',
                        type=int,
                        default=1,
                        choices=range(multiprocessing.cpu_count()),
                        metavar='N',
                        help='number of threads to run (1,{})'.format(multiprocessing.cpu_count()))
    args = parser.parse_args()

    # barcode master list
    print('Loading barcode master list', file=sys.stderr)
    bcm = gsheet.open_by_url('https://docs.google.com/spreadsheets/d/1pxOR4OlODtM8tibT1GX5__W8VzmMBeEV-cFqaUCTVGM/edit#gid=2068628286')
    bcm_df = bcm.worksheet_by_title('Barcode Master List').get_as_df()

    # vector database
    print('Loading vector plate database', file=sys.stderr)
    vec = gsheet.open_by_url('https://docs.google.com/spreadsheets/d/1dEvi81AlWi2tmeo-o2MGTpb3U328gnQB9AVRDxlGMXs/edit#gid=1824873706')
    vec_df = vec.worksheet_by_title('Platemaps').get_as_df()

    # santize inputs
    BAD_BCS = set(['', '#N/A', 'NA'])
    bcm_bc  = [('bcm_1', str(x).upper()) for x in bcm_df['barcode'].unique() if x not in BAD_BCS]
    vec_bc1 = [('vec_1', str(x).upper()) for x in vec_df['bc1'].unique() if x not in BAD_BCS]
    vec_bc2 = [('vec_2', str(x).upper()) for x in vec_df['bc2'].unique() if x not in BAD_BCS]

    # test the reverse compliments too
    bcm_bc_rc   = [('bcm_rc_1', rev_comp(x)) for _,x in bcm_bc]
    vec_bc1_rc  = [('vec_rc_1', rev_comp(x)) for _,x in vec_bc1]
    vec_bc2_rc  = [('vec_rc_2', rev_comp(x)) for _,x in vec_bc2]

    # flatten to enable mapping
    bc_to_test = list(itertools.chain.from_iterable([bcm_bc, vec_bc1, vec_bc2, bcm_bc_rc, vec_bc1_rc, vec_bc2_rc]))

    # collapse rest of columns into tsv
    raw = (x.rstrip().split() for x in args.infile)
    in_clean = (x for x in raw if x[-1] not in BAD_BCS)
    jobs = [(x[-1], args.dist, bc_to_test, '\t'.join(x[:-1])) for x in in_clean]

    print('Checking barcodes', file=sys.stderr)
    well_dict = defaultdict(list)
    with multiprocessing.Pool(args.proc) as pool:
        out_raw = pool.imap_unordered(test_bc, jobs, chunksize=50)
        # recombine split wells
        for tmp_dict in out_raw:
            for key, val in tmp_dict.items():
                well_dict[key].extend(val)

    if args.format == 'tidy':
        for well, bc_list in well_dict.items():
            for bc, dist, bc_db, bc_conflict in bc_list:
                print('\t'.join([well, bc, dist, bc_db, bc_conflict]), file=sys.stdout)
    else:
        for well, bc_list in well_dict.items():
            if all('NA' in x[3] for x in bc_list):
                print(f'{well}\tFALSE', file=sys.stdout)
            else:
                print(f'{well}\tTRUE', file=sys.stdout)
