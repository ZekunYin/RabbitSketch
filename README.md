# libSketch

## Getting Started

A Linux system on a recent x86_64 CPU is required.

### Installing 

Step:

First, we should produce a dynamic depot.

```
cd libsketch

mkdir build

cd build

cmake ..

make
```

We will find exist a file named "libsketch.a" in the build.

### Testing

Then we copy libsketch.a into the catalog named "example" and test.

Step:

```
cd ../examples/

cp ../build/libsketch.a .

make

./minhash 21 ~/workspace/file.fa
```

We will get the value of jaccard and distance.

### TODO
- [ ] Add reverse complement to order minhash
- [ ] Order minhash optimization(profiling)
- [ ] Order minhash optimization(xxHash)
- [ ] Order minhash optimization(sketch and compare)
