TARGET=minhash

CXX=g++
CXXFLAGS=-g -O2 -std=c++17 -fpermissive 

INCLUDE=../src

OBJECTS=test.o libsketch.a

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ `gsl-config --cflags --libs` -I$(INCLUDE) -lz
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

#%.o:%.cpp
test.o:test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz



clean:
	rm *.o $(TARGET)
