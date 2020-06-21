# libSketch

## Getting Started

A Linux system on a recent x86_64 CPU is required.

### Installing (C++ interface) 


```bash
cd libsketch
mkdir build
cd build
cmake -DCXX=ON .. -DCMAKE_INSTALL_PREFIX=.
make
```


### Testing (C++)

```bash
cd ../examples/
#default install dir: ../build/
make 
./minhash genome1.fna genome2.fna
```

We will get the value of jaccard and distance.

### PYTHON bind
**pip install:**
```bash
cd libsketch
pip3 install . --user
```
or
```bash
#pypi available
pip3 install rabbitsketch --user
```
**cmake install**
```bash
cd libsketch
mkdir build
cd build
cmake .. #default with pybind support
make
```
**test using bpython or python**
```bash
cd examples
time python3 test.py #fastx is required
```
### TODO
- [x] Add reverse complement to order minhash
- [x] Order minhash optimization(profiling)
- [x] Order minhash optimization(xxHash)
- [x] Order minhash optimization(sketch and compare)
- [x] remove gsl dependency
- [x] wminhash speed problem
- [ ] reorganize wminhash source
- [ ] remove unrelated minhash code / clean minhash code
- [x] add robin-hood-hashing
- [ ] support multiple types of hashing methods
- [ ] portable to windows and osx?
- [x] remove C++17 requirement (now C++ 14)
- [x] cmake using parameters to compile C++ or python bind
- [x] add setup.py for pip install
- [x] complete pybind.h interface (minhash wmh omh done)
- [x] redesign parameters and all interfaces including minhash wmh omh and hll
- [x] add omhismb2019 license to license file
- [ ] add cpu dispatch
- [ ] pypi publish(dealing with win32 and osx), dealing with requirements
- [ ] pybind std::vector (omh result)
- [ ] revise cmake for make install
- [ ] fix wminhash pybind
- [ ] remove unnecessary installed head files
