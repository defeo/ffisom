
INC_DIR = ../..
BUILD_DIR = ../../build
OBJS = $(wildcard $(BUILD_DIR)/*)
SRC = $(wildcard test*c)
TESTS = $(patsubst %.c, u%, $(SRC))
EXEC = $(patsubst %.c, %, $(SRC))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lflint -lm
.PHONY: clean

clean:
	rm -f $(EXEC)

all:	$(TESTS)

utest_%:	$(OBJS)
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi
