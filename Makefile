# CC = gcc -std=c11
# CFLAGS = -g -O2
# CPP = gcc -E
CXX = g++
# CXXCPP = g++ -E -std=c++11
CPPFLAGS = -I./include
CXXFLAGS = -g -O2 -std=c++14 -mpc80
LDFLAGS = -Xlinker -rpath -Xlinker ./lib
LDLIBS = -L./lib -lmpfr -lgmp -lm -lpthread -lirram -lpng


BIN =  test

all: $(BIN)

# compact: compact.cc
# path: path.cc
test: test.cc
test2: test2.cc


# maintainer-clean: distclean
# distclean: clean
# 	rm -f Makefile

clean:
	rm -f $(BIN)

install: