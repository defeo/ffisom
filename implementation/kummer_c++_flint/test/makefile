INC_DIR = ..
BUILD_DIR = ../build
OBJS = $(wildcard $(BUILD_DIR)/*)
SOURCES = $(wildcard test*.cpp)
TESTS = $(patsubst %.cpp, u%, $(SOURCES))
EXEC = $(patsubst %.cpp, %, $(SOURCES))

GCC = g++
GCCFLAGS = -Wall -O3 -g -I$(INC_DIR)
LIBS = -lflint -lmpfr -lgmp 

.PHONY: clean

all:	$(TESTS)

clean:
	rm -f $(EXEC)

utest_%:	$(OBJS)
	if test -e test_$*.cpp; then\
		$(GCC) -o test_$* test_$*.cpp $(OBJS) $(GCCFLAGS) $(LIBS);\
	fi
