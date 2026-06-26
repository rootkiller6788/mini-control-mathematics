CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wno-misleading-indentation -O2 -g -Iinclude -D_GNU_SOURCE
AR = ar
ARFLAGS = rcs
BLD = build
SRC = src
INC = include
TST = tests
EXM = examples
SRCS = $(SRC)/diffeq_solve.c $(SRC)/ztransform.c $(SRC)/stability.c $(SRC)/statespace.c $(SRC)/convolution.c $(SRC)/applications.c $(SRC)/signal.c $(SRC)/advanced.c
OBJS = $(SRCS:$(SRC)/%.c=$(BLD)/%.o)
LIB = $(BLD)/libdiffeq.a
TEST_SRC = $(TST)/test_diffeq.c
TEST_BIN = $(BLD)/test_diffeq
EXE_FIB = $(BLD)/fibonacci
EXE_FLT = $(BLD)/digital_filter
EXE_POP = $(BLD)/population_model
EXE_PID = $(BLD)/pid_discrete

.PHONY: all lib test examples clean count

all: lib test examples

$(BLD):
	mkdir -p $(BLD)

$(BLD)/%.o: $(SRC)/%.c | $(BLD)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

lib: $(LIB)

$(TEST_BIN): $(TEST_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L$(BLD) -ldiffeq -lm -o $@

test: $(TEST_BIN)
	@echo "=== Running Tests ==="
	$(TEST_BIN)
	@echo "=== All tests passed ==="

$(EXE_FIB): $(EXM)/fibonacci.c $(LIB)
	$(CC) $(CFLAGS) $< -L$(BLD) -ldiffeq -lm -o $@

$(EXE_FLT): $(EXM)/digital_filter.c $(LIB)
	$(CC) $(CFLAGS) $< -L$(BLD) -ldiffeq -lm -o $@

$(EXE_POP): $(EXM)/population_model.c $(LIB)
	$(CC) $(CFLAGS) $< -L$(BLD) -ldiffeq -lm -o $@

$(EXE_PID): $(EXM)/pid_discrete.c $(LIB)
	$(CC) $(CFLAGS) $< -L$(BLD) -ldiffeq -lm -o $@

examples: $(EXE_FIB) $(EXE_FLT) $(EXE_POP) $(EXE_PID)

clean:
	rm -rf $(BLD)

count:
	@echo "include/:" && cat $(INC)/*.h | wc -l
	@echo "src/ (.c only):" && cat $(SRC)/*.c | wc -l
	@echo "Total include/ + src/*.c:" && cat $(INC)/*.h $(SRC)/*.c | wc -l
