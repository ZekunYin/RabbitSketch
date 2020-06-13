TARGET=minhash

CXX=g++
CXXFLAGS=-g -O2 -std=c++17

INCLUDE=../src

OBJECTS=test.o librabbitsketch.a

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ -I$(INCLUDE) -lz
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

#%.o:%.cpp
test.o:test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz



clean:
	rm *.o $(TARGET)
