LIBS=-lm

all: exam exam_dbg exam64 exam64_dbg exam_simd_sse41 exam_alias

exam: main.cpp platform.h rans_byte.h
	g++ -o $@ $< -O3 $(LIBS)

exam_dbg: main.cpp platform.h rans_byte.h
	g++ -o $@ $< -g $(LIBS)

exam64: main64.cpp platform.h rans64.h
	g++ --std=c++20 -o $@ $< -O3 $(LIBS)

exam64_dbg: main64.cpp platform.h rans64.h
	g++ --std=c++20 -o $@ $< -g $(LIBS)

exam_simd_sse41: main_simd.cpp platform.h rans_word_sse41.h
	g++ -o $@ $< -O3 -msse4.1 $(LIBS)

exam_alias: main_alias.cpp platform.h rans_byte.h
	g++ -o $@ $< -O3 $(LIBS)

main_test: main_test.cpp rans_byte_test.h
	g++ --std=c++20 -o $@ $< -O3 $(LIBS) -march=native