INC_DIR = ..
SOURCES = $(wildcard test*.cpp)
TESTS = $(patsubst %.cpp, %, $(SOURCES))

CC = g++
CFLAGS = -Wall -O3 -g -I$(INC_DIR) -L$(INC_DIR)
LIBS = -lkummer -lflint -lmpfr -lgmp 

.PHONY: clean

all: $(TESTS)

clean:
	rm -f $(TESTS)

%: %.cpp
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)
