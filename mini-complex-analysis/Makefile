# Makefile for mini-complex-analysis
# Complex Analysis for Control Theory
# L1-L9: Complete engineering implementation

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/complex_number.c       src/complex_analysis.c       src/complex_integration.c       src/complex_series.c       src/complex_control.c       src/complex_transform.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_complex.c
EX_SRC = examples/ex_nyquist.c          examples/ex_pid_design.c          examples/ex_residue.c

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-complex-analysis: all objects built"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex_nyquist.c $(OBJ) -o examples/ex_nyquist $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_pid_design.c $(OBJ) -o examples/ex_pid_design $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_residue.c $(OBJ) -o examples/ex_residue $(LDFLAGS)

clean:
	rm -f src/*.o tests/test_runner examples/ex_nyquist examples/ex_pid_design examples/ex_residue
