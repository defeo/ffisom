INC_DIR = ..
BUILD_DIR = ../build
OBJS = $(wildcard $(BUILD_DIR)/*)
SOURCES = $(wildcard test*.cpp)
TESTS = $(patsubst %.cpp, u%, $(SOURCES))
EXEC = $(patsubst %.cpp, %, $(SOURCES))

GCC = g++
GCCFLAGS = -Wall -O0 -g -march=native -I$(INC_DIR)
LIBS = -lflint -lgmp -lmpfr

.PHONY: clean

all:	$(TESTS)

clean:
	rm -f $(EXEC)

utest_%:	$(OBJS)
	if test -e test_$*.cpp; then\
		$(GCC) -o test_$* test_$*.cpp $(OBJS) $(GCCFLAGS) $(LIBS);\
	fi