MADAI_GSLPATH = /usr/local
MADAI_CPP = /usr/bin/clang++
#MADAI_CPP = /usr/local/bin/clang-omp++
#MADAI_CPP = /usr/bin/g++
#compiler

MADAI_CFLAGS = -std=c++11 -Wall

INCLUDE = -I ../software/include -I /usr/local/include

LIBRARY = -L ../software/lib -L /usr/local/lib

all :
	make -C ../software
	make bw

bw : bw.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o bw bw.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas
