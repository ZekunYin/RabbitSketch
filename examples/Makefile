TARGET=minhash

CXX=g++
CXXFLAGS=-g -O2 -std=c++14

INCLUDE=../build/include
LIBS=../build/lib
OBJECTS=test.o 

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ -I$(INCLUDE) -lz -lrabbitsketch -L$(LIBS)
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

#%.o:%.cpp
test.o:test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz



clean:
	rm *.o $(TARGET)
