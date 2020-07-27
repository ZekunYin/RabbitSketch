#include <stdint.h>
#include "WangHash.h"

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

/*
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
