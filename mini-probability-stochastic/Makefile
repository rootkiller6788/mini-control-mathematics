# Makefile for mini-probability-stochastic
CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -Iinclude -MMD -MP
AR       = ar rcs
LDFLAGS  = -lm

SRCS = src/probability_core.c src/probability_distributions.c src/probability_stats.c src/probability_estimation.c src/probability_information.c src/stochastic_process.c src/monte_carlo.c src/estimation_filtering.c
OBJS = build/probability_core.o build/probability_distributions.o build/probability_stats.o build/probability_estimation.o build/probability_information.o build/stochastic_process.o build/monte_carlo.o build/estimation_filtering.o
DEPS = $(OBJS:.o=.d)
LIB  = build/libprob.a

.PHONY: all test examples clean

all: $(LIB)

$(LIB): $(OBJS)
	@mkdir -p build
	$(AR) $@ $^
	@echo [OK] Library built

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB)
	$(CC) $(CFLAGS) tests/test_all.c -Lbuild -lprob $(LDFLAGS) -o build/test_all
	./build/test_all && echo [PASS] All tests passed || echo [FAIL] Some tests failed

examples: $(LIB)
	$(CC) $(CFLAGS) examples/ex01_gaussian_mle.c -Lbuild -lprob $(LDFLAGS) -o build/ex01
	$(CC) $(CFLAGS) examples/ex02_markov_chain.c -Lbuild -lprob $(LDFLAGS) -o build/ex02
	$(CC) $(CFLAGS) examples/ex03_kalman_tracking.c -Lbuild -lprob $(LDFLAGS) -o build/ex03
	./build/ex01
	./build/ex02
	./build/ex03

clean:
	rm -rf build

-include $(DEPS)
