CC=g++ 
#CC = icc 

CFLAGS = -fopenmp -std=c++11 -Wall 
#CFLAGS = -qopenmp -std=c++0x -Wall -O3


all: KTRUSS

KTRUSS:
	$(CC) $(CFLAGS) -o ktruss ktruss.cpp 

clean:
	rm -f ./*.exe ./*.o
