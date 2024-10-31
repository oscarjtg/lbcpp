# Compiler
CXX = mpic++

# Compiler flags
CXXFLAGS = -Wall -Wextra -pedantic -Iinclude -I/usr/include/gtest -pthread

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
TEST_DIR = test

# Targets
TARGETS = $(BIN_DIR)/RayleighBenard $(BIN_DIR)/test_equilibria

# Source files for each target
SRC_FILES = $(SRC_DIR)/RayleighBenard.cpp $(wildcard $(SRC_DIR)/lattice/*.cpp) $(wildcard $(SRC_DIR)/macroscopic/*.cpp)
TEST_FILES = $(TEST_DIR)/test_equilibria.cpp

# Object files for each target
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(filter-out $(SRC_DIR)/RayleighBenard.cpp, $(SRC_FILES)))
TEST_OBJ_FILES = $(patsubst $(TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(TEST_FILES))

# Default target
all: $(BIN_DIR)/RayleighBenard $(BIN_DIR)/test_equilibria

# Link object files to create the executables
$(BIN_DIR)/RayleighBenard: $(OBJ_FILES) $(OBJ_DIR)/RayleighBenard.o
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BIN_DIR)/test_equilibria: $(TEST_OBJ_FILES) $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgtest -lgtest_main

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean
