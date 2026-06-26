CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
LDFLAGS = -lm

SRCS = src/matrix.c src/vector.c src/decompositions.c src/solvers.c        src/eigen.c src/control_linalg.c src/matrix_exp.c src/norms.c
OBJS = $(SRCS:.c=.o)
TEST_SRCS = tests/test_basic.c tests/test_decompositions.c tests/test_control.c
TEST_BINS = $(TEST_SRCS:.c=)

.PHONY: all clean test

all: $(OBJS)

test: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "Running $$t..."; ./$$t && echo "PASS: $$t" || echo "FAIL: $$t"; done

tests/test_basic: tests/test_basic.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ tests/test_basic.c $(OBJS) $(LDFLAGS)

tests/test_decompositions: tests/test_decompositions.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ tests/test_decompositions.c $(OBJS) $(LDFLAGS)

tests/test_control: tests/test_control.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ tests/test_control.c $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f src/*.o tests/*.o $(TEST_BINS)
