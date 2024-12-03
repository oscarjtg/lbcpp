# Compiler
CXX = g++

# Compiler warning flags
CXXWARN = -Wall -Wextra -Wpedantic -Werror

# Compiler optimisation flags
CXXOPT = -O3 -march=native

# Paths to header, source, and object files.
INCLUDE_PATH = ../../include
SRC_PATH = ../../src
OBJ_PATH = ../../obj

# Compiler flags
CXXFLAGS = $(CXXWARN) $(CXXOPT) -I$(INCLUDE_PATH)

# Name of target executable file.
TARGET = GaussianHill

# Default target.
all: $(TARGET)

# Source and object files
SRCS = $(shell find $(SRC_PATH) -name '*.cpp')
OBJS = $(patsubst $(SRC_PATH)/%.cpp, $(OBJ_PATH)/%.o, $(SRCS))

# Compile LB source files into object files.
$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile target source file into object file.
$(TARGET).o: $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Link object files to create the executable.
$(TARGET): $(TARGET).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Clean up
clean:
	rm -rf $(OBJS) $(TARGET) $(TARGET).o

# Phony targets
.PHONY: all clean