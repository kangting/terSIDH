.PHONY: main clean

BITS?=128

ifndef UINT_IMPL
UINT_IMPL=uint_custom.c
	ifneq ("$(wildcard p${BITS}/uint1600.s)", "")
		UINT_IMPL=$(wildcard p${BITS}/uint1600.*)
	endif
endif

ifndef FP_IMPL
FP_IMPL=fp.c
	ifneq ("$(wildcard p${BITS}/fp1600.s)", "")
		FP_IMPL=$(wildcard p${BITS}/fp1600.*)
	endif
endif


sources = p${BITS}/constants.c
sources += rng.c
sources += ${UINT_IMPL} ${FP_IMPL}
sources += mont.c
sources += tersidh.c
sources += poly.c
sources += steps.c
sources += fpx.c


includes = $(wildcard *.h p${BITS}/*.h)

main: libtersidh.so
	@cc \
		-I ./ \
		-I p${BITS}/ \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-Wno-unused-but-set-variable \
		-march=native -O3 -funroll-loops -fomit-frame-pointer -flto=6 \
		-Wl,-rpath,. \
		-O2 -g -fstack-protector-strong -fsanitize=address,undefined \
		./libtersidh.so \
		main.c \
		-o main -L./ -ltersidh

libtersidh.so: ${includes} ${sources}
	@cc \
		-shared \
		-fPIC -fvisibility=hidden \
		-g \
		-I ./ \
		-I p${BITS}/ \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-Wno-unused-but-set-variable \
		-march=native -O3 -funroll-loops -fomit-frame-pointer -flto=6 \
		$(sources) \
		-o libtersidh.so


debug: ${includes} ${sources} main.c
	cc \
		-I ./ \
		-I p${BITS}/ \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-g \
		$(sources) \
		main.c \
		-o debug


test: ${includes} $(sources) test.c
	@cc \
    	-I ./p128 \
    	-I ./ \
    	-std=c99 -pedantic \
    	-Wall -Wextra \
    	-march=native -O3 \
		$(sources) \
		test.c \
      	-o test_program


clean:
	@rm -f main debug libtersidh.so test_program *.o p*/*.o

%.o: %.s
	$(CC) -c -fPIC $< -o $@
