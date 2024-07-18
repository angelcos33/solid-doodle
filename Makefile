CC=clang++
INCDIRS=-I.
CFLAGS=-Wall -Wextra -g $(INCDIRS)

CFILES=main.cpp
OBJECTS=main.o

BINARY=main

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CC) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -rf $(BINARY) $(OBJECTS)