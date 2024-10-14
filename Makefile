# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -Iinclude

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
INCLUDE_DIR = include

# Targets
TARGETS = $(BIN_DIR)/RayleighBenard

# Source files for each target
SRC_FILES = $(SRC_DIR)/RayleighBenard.cpp $(wildcard $(SRC_DIR)/lattice/*.cpp) $(wildcard $(SRC_DIR)/macroscopic/*.cpp)

# Object files for each target
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Default target
all: $(TARGETS)

# Link object files to create the executables
$(BIN_DIR)/RayleighBenard: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean
