CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -g
LDFLAGS =

# Graph generator
all: graph_generator

graph_generator: graph_generator.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

graph_generator.o: graph_generator.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f graph_generator graph_generator.o

.PHONY: all clean