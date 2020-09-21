#include <stdint.h>
#include <iostream>
#include <cstring>
#include <immintrin.h>
#include "WangHash.h"

using namespace std;

void print_avx512_epi64(__m512i var)
{
    uint64_t val[8];
    _mm512_store_si512((__m512i *)val, var);
    //printf("The value is: \n");
    for(int i = 0; i < 8; ++i)
        printf("The res is: %llu\n", val[i]);
}

#if defined __AVX512F__ && defined __AVX512BW__
void WangHash_x64_AVX512_x8(const void * seq, const void * seqRev, int length, int kmerSize, void * out)
{
    const int kmerNums    = length - kmerSize + 1;
    const int blockNums   = kmerNums / 8;
    const bool handleTail = (kmerNums % 8)? true : false;
    const uint8_t *data   = (const uint8_t *)seq;
    const uint8_t *dataRev= (const uint8_t *)seqRev;
    uint64_t * tempRes    = (uint64_t *)malloc(sizeof(uint64_t) * kmerNums); //TODO: tempRes can be deleted?
    memset(tempRes, 0, sizeof(uint64_t) * kmerNums);
    const uint8_t mask1 = 0x06;
    uint8_t kmerBuf[kmerSize]; //NOTE: To print the kmer.

    //NOTE: min(seq, seqRev) and charactor to uint
    //TODO: Can charactor_to_uint to be vectorized?
    for(int i = 0; i < kmerNums; ++i)
    {
        bool handleRev = (memcmp((uint8_t *)(data + i), (uint8_t *)(dataRev + length - kmerSize - i), kmerSize) <= 0) ? true : false;
        if(handleRev)
        {
            for(int j = 0; j < kmerSize; ++j)
            {
                kmerBuf[j] = (uint8_t)data[i + j]; //NOTE: To print the kmer.
                uint8_t tempVal = (uint8_t)data[i + j];
                tempVal &= mask1;
                tempVal >>= 1;
                tempRes[i] |= (uint64_t)tempVal;
                tempRes[i] <<= 2;
            }
            //printf("The kmer[%d] is: ", i);
            //printf("The kmer is: ");
            //for(int j = 0; j < kmerSize; ++j)
            //    printf("%c", kmerBuf[j]);
            //printf("\n");
        }
        else
        {
            for(int j = 0; j < kmerSize; ++j)
            {
                kmerBuf[j] = (uint8_t)dataRev[length - kmerSize - i + j]; //NOTE: To print the kmer.
                uint8_t tempVal = (uint8_t)dataRev[length - kmerSize - i + j];
                tempVal &= mask1;
                tempVal >>= 1;
                tempRes[i] |= (uint64_t)tempVal;
                tempRes[i] <<= 2;
            }
            //printf("The kmer[%d] is: ", i);
            //printf("The kmer is: ");
            //for(int j = 0; j < kmerSize; ++j)
            //    printf("%c", kmerBuf[j]);
            //printf("\n");
        }
    }
    //for(int i = 0; i < kmerNums; ++i)
        //printf("tempRes[%d] is: %llu\n", i, tempRes[i]);
        //printf("res is: %llu\n", tempRes[i]);

    __m512i vres;
    __m512i mask2 = _mm512_set1_epi64(0xffffffffffffffff);

    //NOTE: Init the vres and calculate WangHash.
    for(int i = 0; i < blockNums; ++i)
    {
        vres = _mm512_setzero_si512();
        vres = _mm512_load_epi64(tempRes + i * 8);
        //print_avx512_epi64(vres);
        //res = (~res) + (res << 21);
        vres = _mm512_add_epi64(_mm512_xor_epi64(vres,mask2), _mm512_slli_epi64(vres, 21));
        //print_avx512_epi64(vres);
        //res = res ^ (res >> 24);
        vres = _mm512_xor_epi64(vres, _mm512_srli_epi64(vres, 24));
        //print_avx512_epi64(vres);
        //res = (res + (res << 3)) + (res << 8);
        vres = _mm512_add_epi64(_mm512_add_epi64(vres, _mm512_slli_epi64(vres, 3)),
                                                    _mm512_slli_epi64(vres, 8));
        //print_avx512_epi64(vres);
        //res = res ^ (res >> 14);
        vres = _mm512_xor_epi64(vres, _mm512_srli_epi64(vres, 14));
        //print_avx512_epi64(vres);
        //res = (res + (res << 2)) + (res << 4);
        vres = _mm512_add_epi64(_mm512_add_epi64(vres, _mm512_slli_epi64(vres, 2)),
                                                    _mm512_slli_epi64(vres, 4));
        //print_avx512_epi64(vres);
        //res = res ^ (res >> 28);
        vres = _mm512_xor_epi64(vres, _mm512_srli_epi64(vres, 28));
        //print_avx512_epi64(vres);
        //res = res + (res << 31);
        vres = _mm512_add_epi64(vres, _mm512_slli_epi64(vres, 31));
        //print_avx512_epi64(vres);

        //_mm512_storeu_si512(&((uint64_t *)out)[i * 8], vres); //NOTE: SegSIV 
        _mm512_storeu_si512((uint64_t *)out + i * 8, vres); 
    }
    if(handleTail)
        for(int i = blockNums * 8; i < kmerNums; ++i)
        {
            //printf("The tempRes[%d] is: %llu\n", i, tempRes[i]);
            tempRes[i] = (~tempRes[i]) + (tempRes[i] << 21);
            tempRes[i] = tempRes[i] ^ (tempRes[i] >> 24);
            tempRes[i] = (tempRes[i] + (tempRes[i] << 3)) + (tempRes[i] << 8);
            tempRes[i] = tempRes[i] ^ (tempRes[i] >> 14);
            tempRes[i] = (tempRes[i] + (tempRes[i] << 2)) + (tempRes[i] << 4);
            tempRes[i] = tempRes[i] ^ (tempRes[i] >> 28);
            tempRes[i] = tempRes[i] + (tempRes[i] << 31);
            
            //printf("The res is: %llu\n", tempRes[i]);
            ((uint64_t *)out)[i] = tempRes[i];
        }
    //printf("The Hash[%d] is: %llu\n", i, ((uint64_t *)out)[i]);
    free(tempRes);
    return;
}

#endif

void WangHash_x64(const void * kmer, int kmerSize, void * out)
{
    const uint8_t *data = (const uint8_t *)kmer;
    const uint8_t mask = 0x06;
    uint64_t res = 0;
    uint8_t kmerBuf[kmerSize]; //NOTE: To print the kmer.

    for(int i = 0; i < kmerSize; ++i)
    {
        kmerBuf[i] = (uint8_t)data[i]; //NOTE: To print the kmer.
        uint8_t tempVal = (uint8_t)data[i];
        tempVal &= mask;
        tempVal >>= 1;
        res |= (uint64_t)tempVal;
        res <<= 2;
    }
    //printf("The kmer is: ");
    //for(int j = 0; j < kmerSize; ++j)
    //    printf("%c", kmerBuf[j]);
    //printf("\n");
    //printf("res is: %llu\n", res);

    //printf("The res is: %llu\n", res);
    res = (~res) + (res << 21);
    //printf("The res is: %llu\n", res);
    res = res ^ (res >> 24);
    //printf("The res is: %llu\n", res);
    res = (res + (res << 3)) + (res << 8);
    //printf("The res is: %llu\n", res);
    res = res ^ (res >> 14);
    //printf("The res is: %llu\n", res);
    res = (res + (res << 2)) + (res << 4);
    //printf("The res is: %llu\n", res);
    res = res ^ (res >> 28);
    //printf("The res is: %llu\n", res);
    res = res + (res << 31);
    //printf("The res is: %llu\n", res);
    //printf("-------------------------\n");

    ((uint64_t *)out)[0] = res;
    return;
}

/*
void WangHash_x64( const void * key, int len, void * out )
{
    const uint8_t * data = (const uint8_t*)key;
    uint8_t  mask = 0x06;
    uint64_t vkey = 0;

    for(int i = 0; i < len; i++)
    {
        uint8_t meri = (uint8_t)data[i];
        meri &= mask;
        meri >>= 1;
        vkey |= (uint64_t)meri;
        vkey <<= 2;
    }

    vkey = (~vkey) + (vkey << 21);
    vkey = vkey ^ (vkey >> 24);
    vkey = (vkey + (vkey << 3)) + (vkey << 8);
    vkey = vkey ^ (vkey >> 14);
    vkey = (vkey + (vkey << 2)) + (vkey << 4);
    vkey = vkey ^ (vkey >> 28);
    vkey = vkey + (vkey << 31);

    ((uint64_t *)out)[0] = vkey;
}
uint64_t hash_to_uint (const char * kmer, int len）
{
    uint8_t mask = 0x06;
    uint64_t res = 0;
    for(int i = 0; i < len; i++)
    {
        uint8_t meri = (uint8_t)kmer[i];
        meri &= mask;
        meri >>= 1;
        res |= (uint64_t)meri;
        res <<= 2;
    }
    return res;
}
*/
