# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall -Wextra -std=c++23

# Optimization flags for the build target
BUILD_CFLAGS = -O3 

# Directories
SRC_DIR = src
BIN_DIR = bin
TESTS_DIR = tests

# Seed control
ifeq ($(SEED),)
  EXTRA_DEFS :=
else
  EXTRA_DEFS := -DSEED=$(SEED)
endif

# Source files
SRCS = $(filter-out $(SRC_DIR)/tests.cpp, $(wildcard $(SRC_DIR)/*.cpp))
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%.o,$(SRCS))

# Executable and test
EXEC = $(BIN_DIR)/main
MAIN = $(BIN_DIR)/main
TEST_EXEC = $(BIN_DIR)/tests

# Default target
all: $(EXEC)

# Compile with debugging information
debug: BUILD_CFLAGS = -g
debug: $(EXEC)

# Compile target
$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(CFLAGS) $(BUILD_CFLAGS) $(EXTRA_DEFS)

# Compile source files into object files
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $(BUILD_CFLAGS) $(EXTRA_DEFS) -c $< -o $@

# Clean target to remove object files and executable
clean:
	rm -f $(BIN_DIR)/*

# Target for building and running tests
test: $(TEST_EXEC)
	$(TEST_EXEC)

# Compile and link tests/tests.cpp into bin/tests
$(TEST_EXEC): $(TESTS_DIR)/tests.cpp | $(BIN_DIR)
	$(CXX) $(CFLAGS) -g $< -o $@

# Run target to compile and execute the program
run: BUILD_CFLAGS = -g
run: $(EXEC)
	$(EXEC) $(ARGS)