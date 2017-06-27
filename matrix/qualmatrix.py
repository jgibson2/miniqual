import numpy as np
from pysam import AlignedSegment
import pyfaidx

class QualMatrix:
	def __init__(self, chromlengths : dict, reference : pyfaidx.Fasta):
		self.qual_matricies = dict()
		self.cov_matricies = dict()
		self.mismatches = dict()
		self.reference = reference
		for chrom, length in chromlengths.items():
			self.qual_matricies[chrom] = np.zeros(length, dtype=np.float16)
			self.cov_matricies[chrom] = np.zeros(length, dtype=np.int32)
			self.mismatches[chrom] = dict()


	def addSegment(self, segment : AlignedSegment):
		#print(segment.reference_name, segment.reference_start, segment.reference_end)
		refseq = self.reference.get_seq(segment.reference_name, segment.reference_start + 1, segment.reference_end) #1-based start coords
		for index, position in enumerate(segment.get_reference_positions()): #0-based coords
			self.qual_matricies[segment.reference_name][position] = ((self.qual_matricies[segment.reference_name][position] \
				* self.cov_matricies[segment.reference_name][position]) \
				+ segment.query_alignment_qualities[index]) \
				/ (self.cov_matricies[segment.reference_name][position] + 1)
			self.cov_matricies[segment.reference_name][position] += 1
			#print(segment.query_alignment_sequence[index], refseq.seq[position-segment.reference_start])
			if not segment.query_alignment_sequence[index] == refseq.seq[position-segment.reference_start]: #check if the segment base is equal to the reference base
				if position in self.mismatches[segment.reference_name]:
					self.mismatches[segment.reference_name][position] += 1
				else:
					self.mismatches[segment.reference_name][position] = 1
			#print(segment.query_alignment_sequence[index], refseq[position-segment.reference_start])
		#print(refseq)
		#print(segment.query_alignment_sequence)

	def getAverageQuality(self, chrom : str, position : int):
		try:
			return self.qual_matricies[chrom][position]
		except KeyError:
			return 0.0
	
	def getCoverage(self, chrom : str, position : int):
		try:
			return self.cov_matricies[chrom][position]
		except KeyError:
			return 0

	def getQualityMatrix(self, chrom : str):
		return self.qual_matricies[chrom] #intentionally throws KeyError if chrom is not found
	
	def getCoverageMatrix(self, chrom : str):
		return self.cov_matricies[chrom]
