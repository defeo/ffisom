BUILD_DIR = build
SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

CC=g++
CCFLAGS=-Wall -g -O3 -fPIC
LIBS = -lflint -lmpfr -lgmp
LIBRARY=libkummer_nmod.so

.PHONY: clean

all: $(LIBRARY) $(OBJS)

$(LIBRARY): $(OBJS)
	$(CC) -shared $(LDFLAGS) $(LIBS) $(OBJS) -o $@

clean:
	rm -f *.so
	rm -f $(OBJS)
	rm -r $(BUILD_DIR)

$(BUILD_DIR)/%.o:	%.cpp
	echo $(INC_DIR)
	$(CC) -c $< -o $@ $(CCFLAGS)

# make the build dir if not exist
$(OBJS) : | $(BUILD_DIR)
$(BUILD_DIR) :
	mkdir $(BUILD_DIR)
