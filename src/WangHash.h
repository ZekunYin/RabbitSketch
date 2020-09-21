#ifndef _WANGHASH_H_
#define _WANGHASH_H_

void WangHash_x64( const void * key, int len, void * out );
void WangHash_x64_AVX512_x8(const void * seq, const void * seqRev, int length, int kmerSize, void * out);

#endif //_WANGHASH_H_
