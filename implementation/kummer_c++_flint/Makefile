BUILD_DIR = build
SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

CC=g++ -std=c++11
CCFLAGS=-Wall -g -O3 -fPIC
LIBS = -lflint -lmpfr -lgmp
LIBRARY=libkummer.so

.PHONY: clean

all: $(LIBRARY)

$(LIBRARY): $(OBJS)
	$(CC) -shared $(LDFLAGS) $(LIBS) $(OBJS) -o $@

clean:
	rm -f *.so
	rm -f $(OBJS)
	rm -rf $(BUILD_DIR)

$(BUILD_DIR)/%.o:	%.cpp
	$(CC) -c $< -o $@ $(CCFLAGS)

# make the build dir if not exist
$(OBJS) : | $(BUILD_DIR)
$(BUILD_DIR) :
	mkdir $(BUILD_DIR)
