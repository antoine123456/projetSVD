# `make` to build all executables
# `make check` to compile execute tests
# `make progname ; ./progname` to build and run progname.c executable
# `make clean` clean all executables generated

CC=gcc
CFLAGS=-g -Wall
LDFLAGS=-lm -llapacke -fopenmp

SRCS=$(wildcard src/*.c)
OBJS=$(patsubst src/%.c,obj/%.o,$(SRCS))

# Test executables
TEST_SRCS=$(wildcard test/*.c)
TESTS=$(patsubst test/%.c,test/%,$(TEST_SRCS))

# Executables (root-dir)
PROG_SRCS=$(wildcard *.c)
PROGS=$(patsubst %.c,%,$(PROG_SRCS))

# Include directori
INC=-Ihdr

all: $(PROGS) $(TESTS)

# Cancel built-in rule
%: %.c

# Generate programs with `make progname`
%: %.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $^ -o $@ $(LDFLAGS)

# Generate test executables with `make test/testname`
test/%: test/%.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $^ -o $@ $(LDFLAGS)

# Generate every object file
obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ $(LDFLAGS)	

# Clean
clean:
	rm -rf bin/* obj/* $(TESTS) $(PROGS)

# Run every executable tests	
check: $(TESTS)
	@echo
	for test in $(TESTS) ; do ./$$test ; done
	@echo

.PHONY: clean check