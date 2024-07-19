CXX=clang++
INCDIRS=-I.
CXXFLAGS=-Wall -Wextra -g $(INCDIRS) $(DEPFLAGS)
DEPFLAGS=-MP -MD

CFILES=main.cpp
OBJECTS=main.o
DEPFILES=main.d

BINARY=main

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean:
	rm -rf $(BINARY) $(OBJECTS) $(DEPFILES)

-include $(DEPFILES)

.PHONY: all clean