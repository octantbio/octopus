import re
import sys
import argparse
from pathlib import Path

def clean_name(path):
    """
    Assumes names of form path/to/foo_..._R[1,2]...fastq.*
    """
    clean_path = Path(path)
    head = clean_path.name.split('_')[0]
    read = re.search(r"_(R[1-2])_", clean_path.name).group(1)
    ext = ''.join(clean_path.suffixes)
    return head + "_" + read + ext

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create symlinks to all fastq\'s in a sequencing directory. Just specify where the sequencing run is located relative to the current directory and we\'ll handle the rest!')
    parser.add_argument('in_dir',
                        type=str,
                        help='path to sequencing run (or stdin if none)')
    parser.add_argument('-o', '--out-dir',
                        dest='out_dir',
                        type=str,
                        help='where to drop the symlinks (default = current directory)',
                        default=''
                        )
    args = parser.parse_args()

    # dump links in to a folder with the run id
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # grab the reference fasta
    path_to_reference = in_dir.joinpath('input.fasta').resolve()
    if path_to_reference.exists():
        out_dir.joinpath('input.fasta').symlink_to(path_to_reference)
    else:
        print('NO INPUT REFERENCE!\nWill run through de novo assembly only.', file=sys.stdout)

    # grab the SampleSheet
    path_to_samplesheet = in_dir.joinpath('SampleSheet.csv').resolve()
    if path_to_samplesheet.exists():
        out_dir.joinpath('SampleSheet.csv').symlink_to(path_to_samplesheet)
    else:
        raise FileNotFoundError(f'{str(path_to_samplesheet)} does not exist!!')

    # grab the fastqs and drop them into out_dir/fastqs/*.fastq
    out_dir.joinpath('fastqs').mkdir(exist_ok=True)
    fastqs = in_dir.glob('Data/Intensities/BaseCalls/*.fastq*')
    for fq in fastqs:
        out_name = clean_name(fq)
        out_dir.joinpath('fastqs/' + out_name).symlink_to(fq.resolve())
