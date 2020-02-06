MADAI_GSLPATH = /usr/local
MADAI_CPP = /usr/bin/clang++
#MADAI_CPP = /usr/local/bin/clang-omp++
#MADAI_CPP = /usr/bin/g++
#compiler

MADAI_CFLAGS = -Oz -std=c++14

INCLUDE = -I include -I /usr/local/include

LIBRARY = -L lib -L /usr/local/lib

OBJFILES = build/canonical_subs.o build/resonances.o build/parametermap.o build/random.o build/misc.o build/decayproducts.o build/canonical_io.o build/moments.o build/multlist.o build/blastwave.o build/acceptance.o

lib/libcanonical.a : ${OBJFILES}
	ar -ru lib/libcanonical.a ${OBJFILES}

build/decayproducts.o : include/resonances.h src/decayproducts.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/decayproducts.cc -o build/decayproducts.o

build/canonical_subs.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/canonical_subs.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/canonical_subs.cc -o build/canonical_subs.o

build/moments.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/moments.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/moments.cc -o build/moments.o

build/blastwave.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/blastwave.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/blastwave.cc -o build/blastwave.o

build/acceptance.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/acceptance.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/acceptance.cc -o build/acceptance.o

build/multlist.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/multlist.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/multlist.cc -o build/multlist.o
	
build/canonical_io.o : include/canonical.h include/resonances.h include/misc.h include/constants.h include/parametermap.h include/random.h src/canonical_io.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/canonical_io.cc -o build/canonical_io.o
	
build/resonances.o : include/resonances.h include/parametermap.h src/resonances.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/resonances.cc -o build/resonances.o
	
build/parametermap.o : include/resonances.h include/parametermap.h  src/parametermap.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/parametermap.cc -o build/parametermap.o

build/random.o : include/resonances.h include/random.h src/random.cc 
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/random.cc -o build/random.o

build/misc.o : include/constants.h include/misc.h src/misc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/misc.cc -o build/misc.o

clean:
	rm lib/*.a build/*.o