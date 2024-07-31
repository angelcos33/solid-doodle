CXX=g++
INCDIRS=-I.
CXXFLAGS=-Wall -Wextra -std=c++20 -g $(INCDIRS) $(DEPFLAGS)
DEPFLAGS=-MP -MD

CFILES=main.cpp
OBJECTS=main.o
DEPFILES=main.d

# Linker flags
LDFLAGS = 

# Libraries to link
LDLIBS = 

# Executable
BINARY=main

# Cross-platform settings (from https://github.com/KRMisha/Makefile)
ifeq ($(OS),Windows_NT) # OS is a preexisting environment variable on Windows
	OS = windows
else
	UNAME := $(shell uname -s)
	ifeq ($(UNAME),Darwin)
		OS = macos
	else ifeq ($(UNAME),Linux)
		OS = linux
	else
    	$(error OS not supported by this Makefile)
	endif
endif

# OS-specific settings
ifeq ($(OS),windows)
	# Link libgcc and libstdc++ statically on Windows
	LDFLAGS += -static-libgcc -static-libstdc++

	# Windows 32- and 64-bit common settings
	INCDIRS +=
	LDFLAGS +=
	LDLIBS +=

	ifeq ($(win32),1)
		# Windows 32-bit settings
		INCDIRS +=
		LDFLAGS +=
		LDLIBS +=
	else
		# Windows 64-bit settings
		INCDIRS +=
		LDFLAGS +=
		LDLIBS +=
	endif
else ifeq ($(OS),macos)
	# macOS-specific settings
	INCDIRS +=
	LDFLAGS +=
	LDLIBS +=
else ifeq ($(OS),linux)
	# Linux-specific settings
	INCDIRS +=
	LDFLAGS +=
	LDLIBS +=
endif

# Windows-specific default settings
ifeq ($(OS),windows)
	# Add .exe extension to executable
	EXEC := $(BINARY).exe

	ifeq ($(win32),1)
		# Compile for 32-bit
		CXXFLAGS += -m32
	else
		# Compile for 64-bit
		CXXFLAGS += -m64
	endif
endif

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean:
	rm -rf $(BINARY) $(OBJECTS) $(DEPFILES)

-include $(DEPFILES)

.PHONY: all clean