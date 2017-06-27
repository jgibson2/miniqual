import argparse
from matrix import readmatrix
from plotting import *
import re

def __main__():
	parser = argparse.ArgumentParser(description='Generate quality plots for MinION reads')
	parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
	parser.add_argument('-f', '--fastq', help='Fastq file', required=True)
	parser.add_argument('--quality', action='store_true', help='Ouput quality plots')
	parser.add_argument('--coverage', action='store_true', help='Ouput coverage plots')
	parser.add_argument('--length', action='store_true', help='Ouput length plots')
	args = parser.parse_args()

	fq_rec = ['','','','']
	cur_line = 0

	mat = readmatrix.ReadMatrix()

	with open(args.fastq, 'r') as fq:
		for line in fq:
			fq_rec[cur_line % 4] = fq.readline().strip()
			if cur_line % 4 == 3:
				#print(fq_rec)
				try:
					mat.addRecord(fq_rec[1], fq_rec[3]) # add sequecne and quality to matrix
				except ValueError:
					pass
			cur_line += 1

	if args.quality:
		generateReadQualityPlots(mat, args.output_dir)
	if args.coverage:
		generateReadCoveragePlots(mat, args.output_dir)
	if args.length:
		generateReadLengthPlots(mat, args.output_dir)


if __name__ == '__main__':
	__main__()
