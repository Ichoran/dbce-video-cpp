all: dbde_test

dbde_test: makefile \
  dbde_util.h dbde_util.cpp dbde_util.o \
  dbde_util_test.cpp
	g++ -O3 -std=c++14 -march=corei7 dbde_util_test.cpp dbde_util.o -lm -o dbde_test

dbde_util.o: makefile dbde_util.h dbde_util.cpp
	g++ -c -O3 -std=c++14 -march=corei7 dbde_util.cpp
