
INC_DIR = ..
BUILD_DIR = ../build
SOURCES = $(wildcard *.c)
OBJS = $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lm
.PHONY: clean

all:	$(OBJS)

clean:
	rm -f $(OBJS)

$(BUILD_DIR)/%.o:	%.c
	$(CC) -c $< $(CFLAGS) -o $@

