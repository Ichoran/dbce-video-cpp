all: \
  dbde_test

dbde_test: \
  makefile \
  dbde_util.h \
  dbde_util.cpp \
  dbde_util_test.cpp
	g++ -O3 -march=corei7 dbde_util.cpp dbde_util_test.cpp -o dbde_test
