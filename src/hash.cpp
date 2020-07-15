#include "hash.h"

void sketchHashAlgo(const void * key, int len, uint32_t seed, void * out,
                    void (*pHashAlgo)(const void *, int, uint32_t, void *))
{
	pHashAlgo(key, len, seed, out);
	return;
}
