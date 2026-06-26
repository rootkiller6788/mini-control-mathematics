# Makefile for mini-optimization-basics
# Builds the static library, tests, and examples.
#
# Usage: make all | make test | make clean | make lib

CC       = gcc
CFLAGS   = -Wall -Wextra -O2 -std=c11 -fPIC
LDFLAGS  = -lm
AR       = ar

INCDIR   = include
SRCDIR   = src
TSTDIR   = tests
EXADIR   = examples
BLDDIR   = build

SRCS     = $(wildcard $(SRCDIR)/*.c)
OBJS     = $(patsubst $(SRCDIR)/%.c, $(BLDDIR)/%.o, $(SRCS))
LIB      = $(BLDDIR)/liboptimization.a

TESTS    = $(wildcard $(TSTDIR)/*.c)
TESTBINS = $(patsubst $(TSTDIR)/%.c, $(BLDDIR)/test_%, $(TESTS))
EXAMPLES = $(wildcard $(EXADIR)/*.c)
EXABINS  = $(patsubst $(EXADIR)/%.c, $(BLDDIR)/example_%, $(EXAMPLES))

.PHONY: all clean test lib examples

all: lib test examples

$(BLDDIR):
	mkdir -p $(BLDDIR)

$(BLDDIR)/%.o: $(SRCDIR)/%.c | $(BLDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

lib: $(LIB)

$(BLDDIR)/test_%: $(TSTDIR)/%.c $(LIB) | $(BLDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BLDDIR) -loptimization $(LDFLAGS) -o $@

test: $(TESTBINS)
	@echo "=== Running Tests ==="
	@passed=0; failed=0; 	for t in $(TESTBINS); do 		echo "--- $$(basename $$t) ---"; 		if ./$$t; then passed=$$((passed+1)); 		else failed=$$((failed+1)); fi; 	done; 	echo "=== $$passed passed, $$failed failed ==="; 	[ $$failed -eq 0 ] || exit 1

$(BLDDIR)/example_%: $(EXADIR)/%.c $(LIB) | $(BLDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BLDDIR) -loptimization $(LDFLAGS) -o $@

examples: $(EXABINS)

clean:
	rm -rf $(BLDDIR)
