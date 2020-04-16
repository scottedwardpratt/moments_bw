MADAI_GSLPATH = /usr/local
MADAI_CPP = /usr/bin/clang++
#MADAI_CPP = /usr/local/bin/clang-omp++
#MADAI_CPP = /usr/bin/g++
#compiler

MADAI_CFLAGS = -Oz -std=c++11 -Wall

INCLUDE = -I ../software/include -I /usr/local/include

LIBRARY = -L ../software/lib -L /usr/local/lib

all :
	make -C ../software
	make bw
	make cleymans
	make bw_alt
#make canonical

canonical : canonical.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o canonical canonical.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

bw : bw.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o bw bw.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

ideal : ideal.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o ideal ideal.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

bw_alt : bw_alt.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o bw_alt bw_alt.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

cleymans : cleymans.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o cleymans cleymans.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

writemults : writemults.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o writemults writemults.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

analyzemults : analyzemults.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o analyzemults analyzemults.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas

mtester : mtester.cc ../software/lib/libcanonical.a
	${MADAI_CPP} -o mtester mtester.cc ${MADAI_CFLAGS} ${LIBRARY} ${INCLUDE} -lcanonical -lgsl -lgslcblas
