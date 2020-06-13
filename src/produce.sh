swig -c++ -python libsketch.i
gcc -O2 -fPIC -c xxhash.c
g++ -O2 -fPIC -c HashList.cpp HashPriorityQueue.cpp HashSet.cpp MinHashHeap.cpp MurmurHash3.cpp Sketch.cpp countMin.cpp cws.cpp distance.cpp histoSketch.cpp jumpHash.cpp kmerSpectrum.cpp minimizer.cpp 
g++ -O2 -fPIC -c libsketch_wrap.cxx -I/home/old_home/qzh/local/python3/python-3.6.8rc1/include/python3.6m
g++ -shared *.o -o _libsketch.so `gsl-config --cflags --libs` -lz
