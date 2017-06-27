import numpy as np
from pysam import AlignedSegment
import pyfaidx
from math import ceil

class ReadMatrix:
	def __init__(self):
		self.num_reads = 0
		self.max_length = 1
		self.coverage = np.zeros(1, dtype=np.float32)
		self.quality = np.zeros(1, dtype=np.int64)
		self.read_lengths = dict()

	def addRecord(self, sequence : str, quality : str, phred=33):
		if not len(sequence) == len(quality):
			raise ValueError('Sequence length and quality string length are not the same!')
		self.read_lengths[len(sequence)] = self.read_lengths.setdefault(len(sequence), 0) + 1
		if len(sequence) > self.max_length:
			self.quality = np.append(self.quality, [0.0 for i in range(0, ceil((len(sequence)-self.max_length) * 1.5))])
			self.coverage = np.append(self.coverage, [0 for i in range(0, ceil((len(sequence)-self.max_length) * 1.5))])
			self.max_length = self.quality.size
		for index, qual in enumerate(quality):
			self.quality[index] = ((self.quality[index] * self.coverage[index]) + (ord(qual)-phred)) / (self.coverage[index] + 1)
			self.coverage[index] += 1

	def getAverageQuality(self, position : int):
		try:
			return self.quality[position]
		except KeyError:
			return 0.0
	
	def getCoverage(self, position : int):
		try:
			return self.coverage[position]
		except KeyError:
			return 0
