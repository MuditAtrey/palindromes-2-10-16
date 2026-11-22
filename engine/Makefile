# tripal - triple palindrome search
#
#   make            build ./tripal (portable -O3)
#   make native     build with -march=native / -mcpu=apple-m4
#   make test       run the full self-test suite
#   make bench      run the pruning benchmark table
#   make clean

CC      ?= cc
UNAME   := $(shell uname -s)
ARCH    := $(shell uname -m)

CFLAGS  ?= -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter -pthread -Iinclude
LDFLAGS ?= -pthread
LDLIBS  ?= -lm

SRC := src/bigint.c src/orbits.c src/search.c src/checkpoint.c src/main.c
OBJ := $(SRC:.c=.o)
HDR := include/tripal.h include/viable_outer.h

# FLINT (optional): needed only for --dual, the dual-lattice certificate
ifneq ($(wildcard /opt/homebrew/include/flint/fmpz_lll.h),)
DUALF := -DTP_FLINT -I/opt/homebrew/include
DUALL := -L/opt/homebrew/lib -lflint -lgmp -lmpfr
else ifneq ($(wildcard /usr/include/flint/fmpz_lll.h),)
DUALF := -DTP_FLINT
DUALL := -lflint -lgmp -lmpfr
else ifneq ($(wildcard /usr/local/include/flint/fmpz_lll.h),)
DUALF := -DTP_FLINT -I/usr/local/include
DUALL := -L/usr/local/lib -lflint -lgmp -lmpfr
endif

all: tripal

tripal: $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS) $(DUALL)

%.o: %.c $(HDR)
	$(CC) $(CFLAGS) $(DUALF) -c -o $@ $<

# --- optimised build for the host machine -------------------------------
native:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O3 -std=c11 -Wall -Wextra -Wno-unused-parameter -pthread -Iinclude $(NATIVE_FLAGS)" \
	        LDFLAGS="-pthread $(NATIVE_LDFLAGS)" tripal
ifeq ($(UNAME),Darwin)
NATIVE_FLAGS = -mcpu=apple-m4 -funroll-loops
else
NATIVE_FLAGS = -march=native -funroll-loops
endif

# --- tests --------------------------------------------------------------
test: tripal
	python3 tests/run_tests.py

quick: tripal
	python3 tests/run_tests.py --quick

bench: tripal
	@echo
	@echo "=== pruning benchmark (single thread, serial DFS) ==="
	./tripal bench --bits 48,80,128,152,184,192,208,216

clean:
	rm -f $(OBJ) tripal

.PHONY: all native test quick bench clean
