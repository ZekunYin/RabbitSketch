![RabbitSketch](sketch.png)

## Getting Started

A Linux system on a recent x86_64 CPU is required.

### Installing (C++ interface) 


```bash
cd RabbitSketch
mkdir build
cd build
cmake -DCXXAPI=ON .. -DCMAKE_INSTALL_PREFIX=.
make
make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
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
cd RabbitSketch
pip3 install . --user
```
or
```bash
#pypi available (not up to date)
#pip3 install rabbitsketch --user
```
**cmake install**
```bash
cd RabbitSketch
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
