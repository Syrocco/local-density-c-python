 
# Makefile for compiling Pybind11 C++ extension

# Compiler
CXX = g++

# Flags
CXXFLAGS = -O3 -Wall -shared -std=c++11 -fPIC

# Pybind11 includes
PYTHON_INCLUDES := $(shell python3 -m pybind11 --includes)

# Python extension suffix
PYTHON_EXTENSION_SUFFIX := $(shell python3-config --extension-suffix)

# Source files
SRC = cellList.cpp

# Target name
TARGET = cell_list$(PYTHON_EXTENSION_SUFFIX)

# Makefile targets
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) $^ -o $@

clean:
	rm -f $(TARGET)
