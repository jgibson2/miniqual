#ifndef __KMER_H__
#define __KMER_H__

unsigned long kmer_to_num(char* kmer, unsigned int length);
unsigned long sequence_to_num(char* kmer, unsigned int inclusive_start, unsigned int noninclusive_end);
char* num_to_kmer(unsigned long kmer, int length);
unsigned int* count_kmers(char* sequence, unsigned long length, unsigned int k);

#endif
