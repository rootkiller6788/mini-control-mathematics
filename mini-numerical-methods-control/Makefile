CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -O2 -g -Iinclude -c
LDFLAGS = -lm
LINKFLAGS = -lm

SRC_FILES = numerical_core integration linear_systems eigen_methods optimization interpolation control_solvers nonlinear_methods

.PHONY: all clean test examples

all: build/libnumerical.a examples

build/libnumerical.a: $(addprefix src/, $(addsuffix .c, $(SRC_FILES)))
	@mkdir -p build
	for f in $(SRC_FILES); do 		$(CC) $(CFLAGS) src/$$f.c -o build/$$f.o; 	done
	ar rcs build/libnumerical.a build/*.o

test: build/libnumerical.a
	$(CC) -std=c99 -Iinclude tests/test_numerical.c build/libnumerical.a $(LINKFLAGS) -o build/test_numerical
	@echo "Running tests..."
	./build/test_numerical.exe

examples: build/libnumerical.a
	$(CC) -std=c99 -Iinclude examples/example_lqr_controller.c build/libnumerical.a $(LINKFLAGS) -o build/example_lqr
	$(CC) -std=c99 -Iinclude examples/example_ode_integration.c build/libnumerical.a $(LINKFLAGS) -o build/example_ode
	$(CC) -std=c99 -Iinclude examples/example_pole_placement.c build/libnumerical.a $(LINKFLAGS) -o build/example_pole
	@echo "Examples built in build/"

clean:
	rm -rf build/

.DEFAULT_GOAL := all
