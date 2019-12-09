CXX=g++
FLAGS=-O3 -funroll-loops -pipe -march=native -fpermissive -I. -fpic -std=c++17 -Wall -Wextra -Wdisabled-optimization -DNDEBUG -DTHREADSAFE 

ifeq ($(shell uname),Darwin)
	FLAGS := $(FLAGS) -Wa,-q
endif

all: test 

%.o: %.cpp
	$(CXX) -c $(FLAGS)	$< -o $@
%.o: %.c
	$(CXX) -c $(FLAGS)	$< -o $@

test: test.cpp  
	$(CXX) $(FLAGS)	  $< -o $@

clean:
	rm -f test.o test 
