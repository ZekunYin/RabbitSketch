TARGET=minhash

CXX=g++
CXXFLAGS=-g -O2 -std=c++17 -fpermissive 

INCLUDE=../src
#SRC=HashList.cpp HashPriorityQueue.cpp HashSet.cpp MinHashHeap.cpp MurmurHash3.cpp Sketch.cpp test.cpp

OBJECTS=test.o libsketch.a

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ `gsl-config --cflags --libs` -I$(INCLUDE)
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

#%.o:%.cpp
test.o:test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE)



clean:
	rm *.o $(TARGET)
