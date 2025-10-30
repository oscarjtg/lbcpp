# Makefile - build & run C++ tests (Google Test)
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Iinclude -Wall -Wextra
LDFLAGS := -pthread

# Find all prebuilt object files under obj/ (if any)
LIB_OBJS := $(shell find obj -name '*.o' 2>/dev/null)

TEST_SRC := tests/test_equilibria.cpp
TEST_OBJ := $(patsubst %.cpp,%.o,$(TEST_SRC))
BIN_DIR := bin
TEST_BIN := $(BIN_DIR)/test_equilibria

.PHONY: all build test clean gtest-build

all: build

# Build test executable (links library objects + gtest)
build: $(TEST_BIN)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Link step: prefer system gtest libraries; CI should ensure gtest libs are present
$(TEST_BIN): $(BIN_DIR) $(LIB_OBJS) $(TEST_OBJ)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(LIB_OBJS) $(TEST_OBJ) -L/usr/lib -lgtest -lgtest_main $(LDFLAGS)

# Compile test source into object
tests/%.o: tests/%.cpp
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@

# Run tests
test: build
	@echo "Running tests..."
	@./$(TEST_BIN)

# Helper to build system Google Test (Ubuntu installs libgtest-dev as source in /usr/src/gtest)
# Note: building here requires cmake and make; in CI run this target after apt-get install libgtest-dev cmake
gtest-build:
	@if [ -d /usr/src/gtest ]; then \
	  echo "Building Google Test from /usr/src/gtest"; \
	  cd /usr/src/gtest && cmake CMakeLists.txt && make && sudo cp *.a /usr/lib/ || true; \
	else \
	  echo "/usr/src/gtest not found - skip"; \
	fi

clean:
	@echo "Cleaning..."
	-rm -f $(TEST_OBJ) $(TEST_BIN)
	-find . -type f -name '*~' -delete
