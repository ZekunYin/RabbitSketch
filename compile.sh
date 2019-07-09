#!/bin/bash
for cpp in `ls *.cpp`
do 
	g++ -g -c $cpp 
done

g++ -g *.o -o sketch
