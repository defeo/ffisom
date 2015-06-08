BUILD_DIR = build
SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

GCC = g++
GCCFLAGS = -Wall -O0 -g -march=native
LIBS = -lflint -lmpfr -lgmp

.PHONY: clean

all:	$(LIBRARY) $(OBJS)

$(LIBRARY): $(OBJS)
	$(GCC) -shared $(LDFLAGS) $(LIBS) $(OBJS) -o $@

clean:
	rm -f *.so
	rm -f $(OBJS)
	rm -r $(BUILD_DIR)

$(BUILD_DIR)/%.o:	%.cpp
	$(GCC) -c $< -o $@ $(GCCFLAGS)

# make the build dir if not exist
$(OBJS) : | $(BUILD_DIR)
$(BUILD_DIR) :
	mkdir $(BUILD_DIR)
