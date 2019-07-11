TARGET=minhash

CXX=g++
CXXFLAGS=-g 

INCLUDE=../include
#SRC=HashList.cpp HashPriorityQueue.cpp HashSet.cpp MinHashHeap.cpp MurmurHash3.cpp Sketch.cpp test.cpp

OBJECTS=HashList.o HashPriorityQueue.o HashSet.o MinHashHeap.o MurmurHash3.o Sketch.o test.o

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ `gsl-config --cflags --libs`
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



clean:
	rm *.o $(TARGET)
