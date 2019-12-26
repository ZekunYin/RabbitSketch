CXX=g++
#CXX=icpc
#FLAGS=-O3 -funroll-loops -pipe -march=native -fpermissive -I. -fpic -std=c++17 -Wall -Wextra -Wdisabled-optimization -DNDEBUG -DTHREADSAFE 
FLAGS=-O3 -funroll-loops -pipe -march=native -fpermissive -I. -fpic -std=c++17  -Wextra  

ifeq ($(shell uname),Darwin)
	FLAGS := $(FLAGS) -Wa,-q
endif

all: test 

%.o: %.cpp
	$(CXX) -c $(FLAGS)	$< -o $@ -lz
%.o: %.c
	$(CXX) -c $(FLAGS)	$< -o $@ -lz

test: test.cpp  
	$(CXX) $(FLAGS)	  $< -o $@ -lz

clean:
	rm -f test.o test 
