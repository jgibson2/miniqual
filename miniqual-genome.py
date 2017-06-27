import pysam
import pyfaidx
import argparse
from matrix import qualmatrix
from plotting import *
import re

def __main__():
	parser = argparse.ArgumentParser(description='Generate data plots for aligned MinION reads')
	parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
	parser.add_argument('-a', '--alignment', help='Alignment File', required=True)
	parser.add_argument('-t','--type',help="Type of alignment file, 'sam', 'bam', or 'cram'",action='store',default='bam',choices=['bam','sam','cram'])
	parser.add_argument('-r', '--reference', help='Reference FASTA file', required=True)
	parser.add_argument('--regions', nargs='+', help='All regions to output in chr:start-end format with 0-based coordinates, or a chromosome name')
	parser.add_argument('--quality', action='store_true', help='Output quality plots')
	parser.add_argument('--coverage', action='store_true', help='Output coverage plots')
	parser.add_argument('--mismatch', action='store_true', help='Output mismatch plots')
	parser.add_argument('--mapq', action='store', type=int, default=5, help='Minimum mapq score (default 5)')
	args = parser.parse_args()

	alignment_file = None
	if args.type == 'bam':
		alignment_file = pysam.AlignmentFile(args.alignment, 'rb')
	elif args.type == 'cram':
		alignment_file = pysam.AlignmentFile(args.alignment, 'rc')
	else:
		alignment_file = pysam.AlignmentFile(args.alignment, 'r')

	reference_file = pyfaidx.Fasta(args.reference, read_ahead=10000) #much faster when dealing with sequential access (bam/sam sorted by coordinate)

	seq_lengths = dict()
	for chrom in reference_file.keys():
		seq_lengths[chrom] = reference_file.faidx.index[chrom]['rlen']

	mat = qualmatrix.QualMatrix(seq_lengths, reference_file)

	for seg in alignment_file:
		if seg.mapping_quality > args.mapq:
			mat.addSegment(seg)
	if args.regions:
		for region in args.regions:
			if re.match('\w+:\d+\-\d+', region):
				r = re.split('[:\-]', region)
				if not len(r) == 3:
					raise ValueError('Region not formatted correctly')
				if args.quality:
					generateQualityPlots(mat, args.output_dir, chr=r[0], start=int(r[1]), end=int(r[2]))
				if args.coverage:
					generateCoveragePlots(mat, args.output_dir, chr=r[0], start=int(r[1]), end=int(r[2]))
				if args.mismatch:
					generateMismatchPlots(mat, args.output_dir, chr=r[0], start=int(r[1]), end=int(r[2]))
			else:
				if args.quality:
					generateQualityPlots(mat, args.output_dir, chr=region)
				if args.coverage:
					generateCoveragePlots(mat, args.output_dir, chr=region)
				if args.mismatch:
					generateMismatchPlots(mat, args.output_dir, chr=region)
				
	else:
		if args.quality:
			generateQualityPlots(mat, args.output_dir)
		if args.coverage:
			generateCoveragePlots(mat, args.output_dir)
		if args.mismatch:
			generateMismatchPlots(mat, args.output_dir)

if __name__ == '__main__':
	__main__()
