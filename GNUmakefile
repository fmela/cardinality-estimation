CC := $(shell which clang || which gcc)

CXXFLAGS = -Wall -W -O2 -std=c++11 -Wunused -Wno-deprecated-declarations
LDFLAGS = -lm -lstdc++

SRC = $(wildcard *.c *.cc)
HDR = $(wildcard *.h)
BUILD_DIR = bin
BIN = $(SRC:%.c=$(BUILD_DIR)/%) $(SRC:%.cc=$(BUILD_DIR)/%)

all : $(BUILD_DIR) $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%: %.cc $(HDR)
	$(CC) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean
clean:
	if test -d $(BUILD_DIR); then rm -r $(BUILD_DIR); fi
