# ********     Author: Lijun Chang    ******
# ******** Email: 220110422@hit.stu.edu.cn ******
#
CC=g++ -flto -O3

CFLAGS=-c -I. -std=c++11

all: kDC

kDC: .obj/main.o .obj/Graph.o
	${CC} .obj/main.o .obj/Graph.o -o kDC_normal
	rm .obj/*.o

.obj/main.o: main.cpp
	${CC} ${CFLAGS} -o .obj/main.o main.cpp

.obj/Graph.o: Graph.cpp
	${CC} ${CFLAGS} -o .obj/Graph.o Graph.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
