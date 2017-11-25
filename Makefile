### Options
include ./Makefile.options

# Eigen location
INCLUDE=-I$(EIGEN)
# compiler variable
CC = g++
# Flag
CFLAGS = -Wall $(INCLUDE)

all:	gspan

gspan: main.o gspan.o graph.o getopt.o dfs.o cv.o io.o ismin.o misc.o 
	$(CC) $(CFLAGS)  main.o gspan.o graph.o getopt.o dfs.o cv.o io.o ismin.o misc.o -o GspanLarLasso
main.o: main.cpp gspan.h cv.h getopt.h
	$(CC) $(CFLAGS) -c main.cpp
gspan.o: gspan.cpp gspan.h
	$(CC) $(CFLAGS) -c  gspan.cpp
graph.o: graph.cpp gspan.h
	$(CC) $(CFLAGS) -c graph.cpp
getopt.o: getopt.cpp
	$(CC) $(CFLAGS) -c getopt.cpp
dfs.o: dfs.cpp gspan.h
	$(CC) $(CFLAGS) -c dfs.cpp
cv.o: cv.cpp cv.h
	$(CC) $(CFLAGS) -c cv.cpp
io.o: io.cpp gspan.h
	$(CC) $(CFLAGS) -c io.cpp
ismin.o: ismin.cpp gspan.h
	$(CC) $(CFLAGS) -c ismin.cpp
misc.o: misc.cpp gspan.h
	$(CC) $(CFLAGS) -c misc.cpp

clean:
	rm -f *.o