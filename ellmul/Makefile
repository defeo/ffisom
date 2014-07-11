CC=gcc
CFLAGS=-O2 -g
#CFLAGS=-O3 -fomit-frame-pointer

all: ellmul

ellmul: ellmul.c
	${CC} ${CFLAGS} -I${LOCAL}/include -L${LOCAL}/lib $< -lflint -lmpfr -lgmp -lpthread -o $@

clean:
	rm -rf ellmul
