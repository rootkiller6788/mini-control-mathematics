# Makefile — mini-laplace-z-transform
# Compiles library, tests, and examples.

CC       := gcc
CFLAGS   := -std=c99 -Wall -Wextra -O2 -Iinclude
LDFLAGS  := -lm

SRC_DIR  := src
INC_DIR  := include
TST_DIR  := tests
EXM_DIR  := examples
BLD_DIR  := build

# Source files
LIB_SRCS := $(wildcard $(SRC_DIR)/*.c)
LIB_OBJS := $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.o,$(LIB_SRCS))

# Test files
TST_SRCS := $(wildcard $(TST_DIR)/*.c)
TST_BINS := $(patsubst $(TST_DIR)/%.c,$(BLD_DIR)/%,$(TST_SRCS))

# Example files
EXM_SRCS := $(wildcard $(EXM_DIR)/*.c)
EXM_BINS := $(patsubst $(EXM_DIR)/%.c,$(BLD_DIR)/%,$(EXM_SRCS))

.PHONY: all clean test examples lib

all: lib test examples

# Build directory
$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Library objects
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/*.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link test binaries
$(BLD_DIR)/test_%: $(TST_DIR)/test_%.c $(LIB_OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $< $(LIB_OBJS) $(LDFLAGS) -o $@

# Link example binaries
$(BLD_DIR)/example_%: $(EXM_DIR)/example_%.c $(LIB_OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $< $(LIB_OBJS) $(LDFLAGS) -o $@

lib: $(LIB_OBJS)
	@echo "Library built: $(words $(LIB_OBJS)) object files"

test: $(TST_BINS)
	@echo "Running tests..."
	@passed=0; failed=0; \
	for t in $(TST_BINS); do \
		echo "--- $$t ---"; \
		if $$t; then passed=$$((passed+1)); \
		else failed=$$((failed+1)); fi; \
		echo ""; \
	done; \
	echo "=== $$passed passed, $$failed failed ==="

examples: $(EXM_BINS)
	@echo "Running examples..."
	@for e in $(EXM_BINS); do \
		echo "--- $$e ---"; \
		$$e; \
		echo ""; \
	done

clean:
	rm -rf $(BLD_DIR)

# Individual test targets
test-laplace: $(BLD_DIR)/test_laplace
	$(BLD_DIR)/test_laplace

test-stability: $(BLD_DIR)/test_stability
	$(BLD_DIR)/test_stability

# Individual example targets
example-motor: $(BLD_DIR)/example_motor_control
	$(BLD_DIR)/example_motor_control

example-filter: $(BLD_DIR)/example_digital_filter
	$(BLD_DIR)/example_digital_filter

example-response: $(BLD_DIR)/example_system_response
	$(BLD_DIR)/example_system_response
