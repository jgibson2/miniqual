#include "kmer.h"

#include <stdlib.h>
#include <math.h>

unsigned long kmer_to_num(char* kmer, unsigned int length)
{
	unsigned long result = 0;
	for(int i = 0; i < length; i++)
	{
		result <<= 2;
		result |= (long)(kmer[i] & 6) >> 1L;
	}
	return result;
}


unsigned long sequence_to_num(char* kmer, unsigned int inclusive_start, unsigned int noninclusive_end)
{
	unsigned long result = 0;
	for(int i = inclusive_start; i < noninclusive_end; i++)
	{
		result <<= 2L;
		result |= (long)(kmer[i] & 6) >> 1L;
	}
	return result;
}

char* num_to_kmer(unsigned long kmer, int length)
{
	char character_code_replacement[4] = {'A','C','T','G'};
	char* result = (char*)calloc(length+1, sizeof(char));
	for(int i = length - 1; i >= 0; i--)
	{
		result[i] = character_code_replacement[kmer & 3];
		kmer >>= 2L;
	}
	result[length]='\0'; //null-terminate string
	return result;
}

unsigned int* count_kmers(char* sequence, unsigned long length, unsigned int k)
{
	unsigned int* arr = (unsigned int*)calloc(pow(4, k), sizeof(unsigned int));
	unsigned long mer = kmer_to_num(sequence, k);
	for(unsigned long i = 1; i < length - k; i++)
	{
		mer <<= 2L;
		mer |= (long)(sequence[i] & 6L) >> 1L;
		arr[mer] += 1;
	}
	return arr;
}
