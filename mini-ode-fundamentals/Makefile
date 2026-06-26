# Makefile for mini-ode-fundamentals
# Builds the ODE fundamentals library and runs tests.
#
# Usage:
#   make          - build library
#   make test     - build and run tests
#   make examples - build examples
#   make clean    - remove build artifacts

CC       = gcc
CFLAGS   = -std=c99 -D_GNU_SOURCE -Wall -Wextra -Wno-unused-function -Wno-unused-parameter -O2 -g
INCLUDES = -Iinclude
LIBS     = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

# Library sources (C only)
LIB_SRCS = $(SRC_DIR)/ode_ivp.c \
           $(SRC_DIR)/ode_linear.c \
           $(SRC_DIR)/ode_numerical.c \
           $(SRC_DIR)/ode_nonlinear.c \
           $(SRC_DIR)/ode_stability.c \
           $(SRC_DIR)/ode_applications.c

LIB_OBJS = $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.o,$(LIB_SRCS))
LIB_NAME = $(BLD_DIR)/libodefundamentals.a

TEST_SRCS = $(wildcard $(TST_DIR)/*.c)
TEST_BINS = $(patsubst $(TST_DIR)/%.c,$(BLD_DIR)/test_%,$(TEST_SRCS))

EXAMPLE_SRCS = $(wildcard $(EXM_DIR)/*.c)
EXAMPLE_BINS = $(patsubst $(EXM_DIR)/%.c,$(BLD_DIR)/example_%,$(EXAMPLE_SRCS))

.PHONY: all test examples clean

all: $(LIB_NAME)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Compile library objects
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Build static library
$(LIB_NAME): $(LIB_OBJS)
	ar rcs $@ $^

# Build and run tests
test: $(LIB_NAME) $(TEST_BINS)
	@echo "=== Running tests ==="
	@passed=0; failed=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t ---"; \
		if $$t; then \
			echo "  PASS"; passed=$$((passed+1)); \
		else \
			echo "  FAIL"; failed=$$((failed+1)); \
		fi; \
	done; \
	echo "=== $$((passed+failed)) tests: $$passed passed, $$failed failed ==="

# Link test binaries
$(BLD_DIR)/test_%: $(TST_DIR)/%.c $(LIB_NAME) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLD_DIR) -lodefundamentals $(LIBS) -o $@

# Build examples
examples: $(LIB_NAME) $(EXAMPLE_BINS)

$(BLD_DIR)/example_%: $(EXM_DIR)/%.c $(LIB_NAME) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLD_DIR) -lodefundamentals $(LIBS) -o $@

clean:
	rm -rf $(BLD_DIR)
