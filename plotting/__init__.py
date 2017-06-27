import matplotlib
import matplotlib.pyplot as pyplot
from matrix import qualmatrix, readmatrix
import os

def generateQualityPlots(mat : qualmatrix.QualMatrix, out_dir : str, **kwargs):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	if 'chr' in kwargs:
		if 'start' in kwargs and 'end' in kwargs:
			matrix = mat.getQualityMatrix(kwargs['chr'])
			x = range(kwargs['start'], kwargs['end']+1)
			y = matrix[kwargs['start']:kwargs['end']+1]
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='red', markersize=1)
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Average quality')
			pyplot.savefig('%s_%d-%d_quality.png' % (kwargs['chr'], kwargs['start'], kwargs['end']), bbox_inches='tight')
		else:
			matrix = mat.getQualityMatrix(kwargs['chr'])
			x,y = zip(*enumerate(matrix))
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='red', markersize=1)	
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Average quality')
			pyplot.savefig(kwargs['chr'] + '_quality.png', bbox_inches='tight')
	else:
		for chrom, matrix in mat.qual_matricies.items():
			x, y = zip(*enumerate(matrix))
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='red', markersize=1)
			pyplot.title(chrom)
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Average quality')
			pyplot.savefig(chrom + '_quality.png', bbox_inches='tight')


def generateCoveragePlots(mat : qualmatrix.QualMatrix, out_dir : str, **kwargs):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	if 'chr' in kwargs:
		if 'start' in kwargs and 'end' in kwargs:
			matrix = mat.getCoverageMatrix(kwargs['chr'])
			x = range(kwargs['start'], kwargs['end']+1)
			y = matrix[kwargs['start']:kwargs['end']+1]
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='blue', markersize=1)
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Coverage')
			pyplot.savefig('%s_%d-%d_coverage.png' % (kwargs['chr'], kwargs['start'], kwargs['end']), bbox_inches='tight')
		else:
			matrix = mat.getCoverageMatrix(kwargs['chr'])
			x,y = zip(*enumerate(matrix))
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='blue', markersize=1)
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Coverage')
			pyplot.savefig(kwargs['chr'] + '_coverage.png', bbox_inches='tight')
	else:
		for chrom, matrix in mat.cov_matricies.items():
			x, y = zip(*enumerate(matrix))
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='blue', markersize=1)
			pyplot.title(chrom)
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Coverage')
			pyplot.savefig(chrom + '_coverage.png', bbox_inches='tight')


def generateMismatchPlots(mat : qualmatrix.QualMatrix, out_dir : str, **kwargs):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	if 'chr' in kwargs:
		if 'start' in kwargs and 'end' in kwargs:
			x = [x for x in range(kwargs['start'], kwargs['end']+1)]
			y = [mat.mismatches[kwargs['chr']][i]/mat.cov_matricies[kwargs['chr']][i] if i in mat.mismatches[kwargs['chr']] else 0 for i in x]
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='green', markersize=1)
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Mismatch frequency')
			pyplot.savefig('%s_%d-%d_mismatch.png' % (kwargs['chr'], kwargs['start'], kwargs['end']), bbox_inches='tight')
		else:
			matrix = mat.getCoverageMatrix(kwargs['chr'])
			x = [x for x in range(matrix.size)]
			y = [mat.mismatches[kwargs['chr']][i]/mat.cov_matricies[kwargs['chr']][i] if i in mat.mismatches[kwargs['chr']] else 0 for i in x]
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='green', markersize=1)
			pyplot.title(kwargs['chr'])
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Mismatch frequency')
			pyplot.savefig(kwargs['chr'] + '_mismatch.png', bbox_inches='tight')
	else:
		for chrom, matrix in mat.qual_matricies.items():
			x = [x for x in range(matrix.size)]
			y = [mat.mismatches[chrom][i]/mat.cov_matricies[chrom][i] if i in mat.mismatches[chrom] else 0 for i in x]
			pyplot.figure(figsize=(20,10))
			pyplot.plot(x, y, '-o', color='green', markersize=1)
			pyplot.title(chrom)
			pyplot.xlabel('Genomic coordinate')
			pyplot.ylabel('Mismatch frequency')
			pyplot.savefig(chrom + '_mismatch.png', bbox_inches='tight')


def generateReadQualityPlots(mat : readmatrix.ReadMatrix, out_dir : str):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	x = [i for i in range(mat.max_length)]
	y = [mat.getAverageQuality(i) for i in x]
	pyplot.figure(figsize=(20,10))
	pyplot.plot(x, y, '-o', color='red', markersize=1)
	pyplot.title('Average Qualities')
	pyplot.xlabel('Read coordinate')
	pyplot.ylabel('Average quality')
	pyplot.savefig('read_quality.png', bbox_inches='tight')


def generateReadCoveragePlots(mat : readmatrix.ReadMatrix, out_dir : str):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	x = [i for i in range(mat.max_length)]
	y = [mat.getCoverage(i) for i in x]
	pyplot.figure(figsize=(20,10))
	pyplot.plot(x, y, '-o', color='blue', markersize=1)
	pyplot.title('Read Coverage')
	pyplot.xlabel('Read coordinate')
	pyplot.ylabel('Coverage')
	pyplot.savefig('read_coverage.png', bbox_inches='tight')

def generateReadLengthPlots(mat : readmatrix.ReadMatrix, out_dir : str):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(out_dir)

	x = [i for i in range(mat.max_length + 1)]
	y = [mat.read_lengths.setdefault(i, None) for i in x]
	pyplot.figure(figsize=(20,10))
	pyplot.scatter(x, y, color='green')
	pyplot.title('Read Lengths')
	pyplot.xlabel('Read Length')
	pyplot.ylabel('Number of Reads')
	pyplot.savefig('read_lengths.png', bbox_inches='tight')
