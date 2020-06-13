%module libsketch

%{
#include "HashList.h"
#include "HashPriorityQueue.h"
#include "HashSet.h"
#include "MinHashHeap.h"
#include "MurmurHash3.h"
#include "Sketch.h"
#include "bloom_filter.hpp"
#include "countMin.h"
#include "cws.h"
#include "distance.h"
#include "hash.h"
#include "hash_int.h"
#include "histoSketch.h"
#include "jumpHash.h"
#include "kmerSpectrum.h"
#include "minimizer.h"
#include "xxh3.h"
#include "xxhash.h"
#include "xxhash.hpp"
%}

%include <stdint.i>
%include "Sketch.h"
