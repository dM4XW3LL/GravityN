# ===========================================================================
#  Makefile — GravityN  (N-body Gravitational Simulator)
# ===========================================================================
#
#  Targets:
#    make / make all   Build the shared physics library (nbody.so/.dylib/.dll)
#    make test         Build and run the C validation suite (tests/test_nbody)
#    make run          Build the library (if needed) and launch the GUI
#    make clean        Remove all build artifacts
#    make rebuild      clean, then build everything again
#    make help         Show this help message
#
#  On Windows, run these targets from an MSYS2 / Git Bash / WSL shell so
#  that the Unix-style commands below (rm, PYTHONPATH=...) work correctly.
# ===========================================================================

CC     := gcc
CFLAGS := -O2 -Wall -Wextra
LDLIBS := -lm

SRC_DIR  := src/core
TEST_DIR := tests
GUI_DIR  := src/gui

CORE_SRCS := $(SRC_DIR)/nbody.c $(SRC_DIR)/kepler.c
CORE_HDRS := $(SRC_DIR)/nbody.h $(SRC_DIR)/kepler.h

# ── Platform detection ───────────────────────────────────────────────────
# Selects the shared-library name/flags, executable suffix, and python
# command for the current platform.
ifeq ($(OS),Windows_NT)
    LIB_NAME  := nbody.dll
    LIB_FLAGS := -shared
    LDLIBS    :=
    EXE_EXT   := .exe
    PYTHON    := python
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LIB_NAME  := nbody.dylib
        LIB_FLAGS := -dynamiclib -fPIC
    else
        LIB_NAME  := nbody.so
        LIB_FLAGS := -shared -fPIC
    endif
    EXE_EXT :=
    PYTHON  := python3
endif

TEST_BIN := $(TEST_DIR)/test_nbody$(EXE_EXT)

.PHONY: all lib test run clean rebuild help

# ── Default target ──────────────────────────────────────────────────────
all: lib

# ── Shared physics library ──────────────────────────────────────────────
lib: $(LIB_NAME)

$(LIB_NAME): $(CORE_SRCS) $(CORE_HDRS)
	$(CC) $(CFLAGS) $(LIB_FLAGS) -o $@ $(CORE_SRCS) $(LDLIBS)

# ── C validation suite ───────────────────────────────────────────────────
# tests/test_nbody.c includes "../src/core/nbody.h" directly, so no -I flag
# is needed as long as this is built from the project root.
test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_DIR)/test_nbody.c $(CORE_SRCS) $(CORE_HDRS)
	$(CC) $(CFLAGS) -o $@ $(TEST_DIR)/test_nbody.c $(CORE_SRCS) $(LDLIBS)

# ── GUI ───────────────────────────────────────────────────────────────────
run: lib
	PYTHONPATH=$(GUI_DIR) $(PYTHON) $(GUI_DIR)/nbodysim.py

# ── Housekeeping ─────────────────────────────────────────────────────────
clean:
	rm -f nbody.so nbody.dylib nbody.dll $(TEST_BIN)

rebuild: clean all

help:
	@echo "GravityN -- available make targets:"
	@echo ""
	@echo "  make          Build the shared physics library ($(LIB_NAME))"
	@echo "  make lib      Same as above"
	@echo "  make test     Build and run the C validation suite"
	@echo "  make run      Build the library and launch the GUI"
	@echo "  make clean    Remove all build artifacts"
	@echo "  make rebuild  clean, then build everything again"
