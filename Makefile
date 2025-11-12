# ============================================================================
# FDPricing Library - Makefile
# ============================================================================

# ============================================================================
# Build Configuration
# ============================================================================

# Compiler selection: gcc or clang
CC ?= gcc

# Build mode: debug or release
BUILD_MODE ?= release

# Installation prefix
PREFIX ?= /usr/local

# ============================================================================
# Compiler flags
# ============================================================================

# Common flags
COMMON_FLAGS = -Wall -Wextra -Wpedantic -std=c99 -fPIC
INCLUDES = -Iinclude

# Release mode flags
RELEASE_FLAGS = -O3 -DNDEBUG -march=native
ifeq ($(CC),gcc)
    RELEASE_FLAGS += -flto
endif

# Debug mode flags
DEBUG_FLAGS = -O0 -g -DDEBUG -fsanitize=address -fsanitize=undefined

# Set flags based on build mode
ifeq ($(BUILD_MODE),debug)
    CFLAGS = $(COMMON_FLAGS) $(DEBUG_FLAGS)
    LDFLAGS_EXTRA = -fsanitize=address -fsanitize=undefined
else
    CFLAGS = $(COMMON_FLAGS) $(RELEASE_FLAGS)
    LDFLAGS_EXTRA = 
endif

LDFLAGS = -lm $(LDFLAGS_EXTRA)

# Test-specific flags (add Criterion library)
TEST_LDFLAGS = $(LDFLAGS) -lcriterion

# ============================================================================
# Directories
# ============================================================================

BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
LIB_DIR = $(BUILD_DIR)/lib
BIN_DIR = $(BUILD_DIR)/bin
TEST_DIR = tests
TEST_BIN_DIR = $(BUILD_DIR)/tests

# Object subdirectories
OBJ_CORE_DIR = $(OBJ_DIR)/core
OBJ_MODELS_DIR = $(OBJ_DIR)/models
OBJ_OPTIONS_DIR = $(OBJ_DIR)/options
OBJ_NUMERICS_DIR = $(OBJ_DIR)/numerics
OBJ_UTILS_DIR = $(OBJ_DIR)/utils

# ============================================================================
# Colors for pretty output
# ============================================================================

COLOR_RESET = \033[0m
COLOR_BOLD = \033[1m
COLOR_RED = \033[31m
COLOR_GREEN = \033[32m
COLOR_YELLOW = \033[33m
COLOR_BLUE = \033[34m
COLOR_CYAN = \033[36m

# ============================================================================
# Source files
# ============================================================================

# Core sources
CORE_SRCS = $(wildcard src/core/*.c)
CORE_OBJS = $(patsubst src/core/%.c,$(OBJ_CORE_DIR)/%.o,$(CORE_SRCS))

# Model sources
MODELS_SRCS = $(wildcard src/models/*.c)
MODELS_OBJS = $(patsubst src/models/%.c,$(OBJ_MODELS_DIR)/%.o,$(MODELS_SRCS))

# Options sources
OPTIONS_SRCS = $(wildcard src/options/*.c)
OPTIONS_OBJS = $(patsubst src/options/%.c,$(OBJ_OPTIONS_DIR)/%.o,$(OPTIONS_SRCS))

# Numerics sources
NUMERICS_SRCS = $(wildcard src/numerics/*.c)
NUMERICS_OBJS = $(patsubst src/numerics/%.c,$(OBJ_NUMERICS_DIR)/%.o,$(NUMERICS_SRCS))

# Utils sources
UTILS_SRCS = $(wildcard src/utils/*.c)
UTILS_OBJS = $(patsubst src/utils/%.c,$(OBJ_UTILS_DIR)/%.o,$(UTILS_SRCS))

# All objects
ALL_OBJS = $(CORE_OBJS) $(MODELS_OBJS) $(OPTIONS_OBJS) $(NUMERICS_OBJS) $(UTILS_OBJS)

# Test sources
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS = $(patsubst $(TEST_DIR)/%.c,$(TEST_BIN_DIR)/%,$(TEST_SRCS))

# ============================================================================
# Library targets
# ============================================================================

LIB_NAME = fdpricing
STATIC_LIB = $(LIB_DIR)/lib$(LIB_NAME).a
SHARED_LIB = $(LIB_DIR)/lib$(LIB_NAME).so

# ============================================================================
# Example programs
# ============================================================================

EXAMPLES = simple_european american_comparison grid_convergence
EXAMPLE_BINS = $(addprefix $(BIN_DIR)/,$(EXAMPLES))

# ============================================================================
# Main targets
# ============================================================================

.PHONY: all build lib static shared examples tests clean distclean install uninstall help
.DEFAULT_GOAL := help

# Help target
help:
	@echo "$(COLOR_BOLD)$(COLOR_CYAN)FDPricing Library Build System$(COLOR_RESET)"
	@echo "$(COLOR_BOLD)==============================$(COLOR_RESET)"
	@echo ""
	@echo "$(COLOR_BOLD)Configuration:$(COLOR_RESET)"
	@echo "  CC=$(COLOR_GREEN)$(CC)$(COLOR_RESET)"
	@echo "  BUILD_MODE=$(COLOR_GREEN)$(BUILD_MODE)$(COLOR_RESET)"
	@echo ""
	@echo "$(COLOR_BOLD)Build Targets:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)make build$(COLOR_RESET)         - Build libraries and examples ($(BUILD_MODE) mode)"
	@echo "  $(COLOR_YELLOW)make lib$(COLOR_RESET)           - Build static and shared libraries"
	@echo "  $(COLOR_YELLOW)make static$(COLOR_RESET)        - Build static library only"
	@echo "  $(COLOR_YELLOW)make shared$(COLOR_RESET)        - Build shared library only"
	@echo "  $(COLOR_YELLOW)make examples$(COLOR_RESET)      - Build example programs"
	@echo "  $(COLOR_YELLOW)make tests$(COLOR_RESET)         - Build test suite"
	@echo ""
	@echo "$(COLOR_BOLD)Test Targets:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)make test$(COLOR_RESET)          - Run all tests with Criterion"
	@echo "  $(COLOR_YELLOW)make test-verbose$(COLOR_RESET)  - Run tests with verbose output"
	@echo "  $(COLOR_YELLOW)make test-list$(COLOR_RESET)     - List all available tests"
	@echo "  $(COLOR_YELLOW)make test-quick$(COLOR_RESET)    - Run quick smoke tests"
	@echo ""
	@echo "$(COLOR_BOLD)Run Targets:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)make run-simple$(COLOR_RESET)    - Run simple_european example"
	@echo "  $(COLOR_YELLOW)make run-american$(COLOR_RESET)  - Run american_comparison example"
	@echo "  $(COLOR_YELLOW)make run-convergence$(COLOR_RESET) - Run grid_convergence example"
	@echo "  $(COLOR_YELLOW)make run-all$(COLOR_RESET)       - Run all examples"
	@echo ""
	@echo "$(COLOR_BOLD)Clean Targets:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)make clean$(COLOR_RESET)         - Remove build artifacts"
	@echo "  $(COLOR_YELLOW)make distclean$(COLOR_RESET)     - Remove all generated files"
	@echo ""
	@echo "$(COLOR_BOLD)Install Targets:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)make install$(COLOR_RESET)       - Install library to $(PREFIX)"
	@echo "  $(COLOR_YELLOW)make uninstall$(COLOR_RESET)     - Remove installed library"
	@echo ""
	@echo "$(COLOR_BOLD)Configuration Options:$(COLOR_RESET)"
	@echo "  $(COLOR_YELLOW)CC=clang$(COLOR_RESET)           - Use clang compiler"
	@echo "  $(COLOR_YELLOW)CC=gcc$(COLOR_RESET)             - Use gcc compiler (default)"
	@echo "  $(COLOR_YELLOW)BUILD_MODE=debug$(COLOR_RESET)   - Build with debug symbols and sanitizers"
	@echo "  $(COLOR_YELLOW)BUILD_MODE=release$(COLOR_RESET) - Build with optimizations (default)"
	@echo "  $(COLOR_YELLOW)PREFIX=/path$(COLOR_RESET)       - Set installation prefix"
	@echo ""
	@echo "$(COLOR_BOLD)Example Usage:$(COLOR_RESET)"
	@echo "  make build                      # Build everything (release)"
	@echo "  make build BUILD_MODE=debug     # Build everything (debug)"
	@echo "  make CC=clang build             # Build with clang"
	@echo "  make test                       # Run tests"
	@echo "  make clean build                # Clean rebuild"
	@echo "  sudo make install               # Install system-wide"

# Default 'all' target - just build
all: build

# Main build target
build: info dirs $(STATIC_LIB) $(SHARED_LIB) examples
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_GREEN)✓ Build complete!$(COLOR_RESET)"
	@echo "  Static library:  $(COLOR_CYAN)$(STATIC_LIB)$(COLOR_RESET)"
	@echo "  Shared library:  $(COLOR_CYAN)$(SHARED_LIB)$(COLOR_RESET)"
	@echo "  Examples:        $(COLOR_CYAN)$(BIN_DIR)/$(COLOR_RESET)"
	@echo ""
	@echo "Run '$(COLOR_YELLOW)make test$(COLOR_RESET)' to run tests"
	@echo "Run '$(COLOR_YELLOW)make run-all$(COLOR_RESET)' to run examples"

# Print build configuration
info:
	@echo "$(COLOR_BOLD)$(COLOR_BLUE)Building FDPricing Library$(COLOR_RESET)"
	@echo "  Compiler:   $(COLOR_GREEN)$(CC)$(COLOR_RESET)"
	@echo "  Mode:       $(COLOR_GREEN)$(BUILD_MODE)$(COLOR_RESET)"
	@echo "  CFLAGS:     $(CFLAGS)"
	@echo ""

# Create directories
dirs:
	@mkdir -p $(OBJ_CORE_DIR)
	@mkdir -p $(OBJ_MODELS_DIR)
	@mkdir -p $(OBJ_OPTIONS_DIR)
	@mkdir -p $(OBJ_NUMERICS_DIR)
	@mkdir -p $(OBJ_UTILS_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(TEST_BIN_DIR)

# ============================================================================
# Compilation rules
# ============================================================================

$(OBJ_CORE_DIR)/%.o: src/core/%.c
	@echo "  $(COLOR_CYAN)CC$(COLOR_RESET)      $<"
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_MODELS_DIR)/%.o: src/models/%.c
	@echo "  $(COLOR_CYAN)CC$(COLOR_RESET)      $<"
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_OPTIONS_DIR)/%.o: src/options/%.c
	@echo "  $(COLOR_CYAN)CC$(COLOR_RESET)      $<"
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_NUMERICS_DIR)/%.o: src/numerics/%.c
	@echo "  $(COLOR_CYAN)CC$(COLOR_RESET)      $<"
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_UTILS_DIR)/%.o: src/utils/%.c
	@echo "  $(COLOR_CYAN)CC$(COLOR_RESET)      $<"
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ============================================================================
# Library targets
# ============================================================================

lib: $(STATIC_LIB) $(SHARED_LIB)

static: $(STATIC_LIB)

shared: $(SHARED_LIB)

$(STATIC_LIB): $(ALL_OBJS)
	@echo "  $(COLOR_GREEN)AR$(COLOR_RESET)      $@"
	@ar rcs $@ $(ALL_OBJS)

$(SHARED_LIB): $(ALL_OBJS)
	@echo "  $(COLOR_GREEN)LD$(COLOR_RESET)      $@"
	@$(CC) -shared -o $@ $(ALL_OBJS) $(LDFLAGS)

# ============================================================================
# Example programs
# ============================================================================

examples: $(EXAMPLE_BINS)

$(BIN_DIR)/%: examples/%.c $(STATIC_LIB)
	@echo "  $(COLOR_GREEN)CCLD$(COLOR_RESET)    $@"
	@$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ $(STATIC_LIB) $(LDFLAGS)

# ============================================================================
# Test suite with Criterion
# ============================================================================

tests: $(TEST_BINS)

$(TEST_BIN_DIR)/%: $(TEST_DIR)/%.c $(STATIC_LIB)
	@echo "  $(COLOR_GREEN)CCLD$(COLOR_RESET)    $@ (Criterion)"
	@$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ $(STATIC_LIB) $(TEST_LDFLAGS)

# Check if Criterion is installed
check-criterion:
	@if ! command -v pkg-config >/dev/null 2>&1 || ! pkg-config --exists criterion 2>/dev/null; then \
		if ! ldconfig -p 2>/dev/null | grep -q libcriterion; then \
			echo "$(COLOR_BOLD)$(COLOR_RED)Error: Criterion library not found!$(COLOR_RESET)"; \
			echo ""; \
			echo "Please install Criterion:"; \
			echo "  Ubuntu/Debian: sudo apt-get install libcriterion-dev"; \
			echo "  macOS:         brew install criterion"; \
			echo "  From source:   https://github.com/Snaipe/Criterion"; \
			echo ""; \
			exit 1; \
		fi \
	fi

# Run all tests with Criterion
test: check-criterion tests
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_BLUE)Running Test Suite (Criterion)$(COLOR_RESET)"
	@echo "$(COLOR_BOLD)===============================$(COLOR_RESET)"
	@echo ""
	@passed=0; failed=0; total_tests=0; \
	for test_bin in $(TEST_BINS); do \
		echo "$(COLOR_CYAN)Running: $$(basename $$test_bin)$(COLOR_RESET)"; \
		if $$test_bin --verbose=0; then \
			echo "  $(COLOR_GREEN)✓ PASSED$(COLOR_RESET)"; \
			passed=$$((passed + 1)); \
		else \
			echo "  $(COLOR_RED)✗ FAILED$(COLOR_RESET)"; \
			failed=$$((failed + 1)); \
		fi; \
		echo ""; \
	done; \
	echo "$(COLOR_BOLD)========================================$(COLOR_RESET)"; \
	echo "$(COLOR_BOLD)Test Summary$(COLOR_RESET)"; \
	echo "$(COLOR_BOLD)========================================$(COLOR_RESET)"; \
	echo "Test Suites: $$((passed + failed))"; \
	echo "Passed:      $(COLOR_GREEN)$$passed$(COLOR_RESET)"; \
	echo "Failed:      $(COLOR_RED)$$failed$(COLOR_RESET)"; \
	echo "$(COLOR_BOLD)========================================$(COLOR_RESET)"; \
	if [ $$failed -eq 0 ]; then \
		echo ""; \
		echo "$(COLOR_BOLD)$(COLOR_GREEN)✓ All test suites passed!$(COLOR_RESET)"; \
		echo ""; \
		exit 0; \
	else \
		echo ""; \
		echo "$(COLOR_BOLD)$(COLOR_RED)✗ Some test suites failed!$(COLOR_RESET)"; \
		echo ""; \
		exit 1; \
	fi

# Run tests with verbose Criterion output
test-verbose: check-criterion tests
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_BLUE)Running Test Suite (Verbose)$(COLOR_RESET)"
	@echo "$(COLOR_BOLD)=============================$(COLOR_RESET)"
	@echo ""
	@for test_bin in $(TEST_BINS); do \
		echo "$(COLOR_CYAN)Running: $$(basename $$test_bin)$(COLOR_RESET)"; \
		echo "$(COLOR_BOLD)----------------------------------------$(COLOR_RESET)"; \
		$$test_bin --verbose=1 || exit 1; \
		echo ""; \
	done

# List all tests without running them
test-list: check-criterion tests
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_BLUE)Available Tests$(COLOR_RESET)"
	@echo "$(COLOR_BOLD)===============$(COLOR_RESET)"
	@echo ""
	@for test_bin in $(TEST_BINS); do \
		echo "$(COLOR_CYAN)$$(basename $$test_bin):$(COLOR_RESET)"; \
		$$test_bin --list || true; \
		echo ""; \
	done

# Run a specific test suite
test-suite-%: check-criterion $(TEST_BIN_DIR)/test_%
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_CYAN)Running: test_$*$(COLOR_RESET)"
	@echo ""
	@$(TEST_BIN_DIR)/test_$*

# Quick smoke test (just run test_basic)
test-quick: check-criterion $(TEST_BIN_DIR)/test_basic
	@echo ""
	@echo "$(COLOR_BOLD)$(COLOR_BLUE)Running Quick Tests$(COLOR_RESET)"
	@echo "$(COLOR_BOLD)===================$(COLOR_RESET)"
	@echo ""
	@$(TEST_BIN_DIR)/test_basic

# Run tests with TAP output (for CI systems)
test-tap: check-criterion tests
	@for test_bin in $(TEST_BINS); do \
		$$test_bin --tap; \
	done

# Run tests with JUnit XML output (for CI systems)
test-junit: check-criterion tests
	@mkdir -p $(BUILD_DIR)/test-results
	@for test_bin in $(TEST_BINS); do \
		$$test_bin --xml=$(BUILD_DIR)/test-results/$$(basename $$test_bin).xml; \
	done

# ============================================================================
# Run examples
# ============================================================================

.PHONY: run-simple run-american run-convergence run-all

run-simple: $(BIN_DIR)/simple_european
	@echo "$(COLOR_BOLD)$(COLOR_CYAN)Running: simple_european$(COLOR_RESET)"
	@echo ""
	@$(BIN_DIR)/simple_european

run-american: $(BIN_DIR)/american_comparison
	@echo "$(COLOR_BOLD)$(COLOR_CYAN)Running: american_comparison$(COLOR_RESET)"
	@echo ""
	@$(BIN_DIR)/american_comparison

run-convergence: $(BIN_DIR)/grid_convergence
	@echo "$(COLOR_BOLD)$(COLOR_CYAN)Running: grid_convergence$(COLOR_RESET)"
	@echo ""
	@$(BIN_DIR)/grid_convergence

run-all: run-simple run-american run-convergence

# ============================================================================
# Installation
# ============================================================================

INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_INC_DIR = $(PREFIX)/include

install: build
	@echo "$(COLOR_BOLD)Installing fdpricing library...$(COLOR_RESET)"
	@mkdir -p $(INSTALL_LIB_DIR)
	@mkdir -p $(INSTALL_INC_DIR)
	@cp $(STATIC_LIB) $(INSTALL_LIB_DIR)/
	@cp $(SHARED_LIB) $(INSTALL_LIB_DIR)/
	@cp include/fdpricing.h $(INSTALL_INC_DIR)/
	@ldconfig 2>/dev/null || true
	@echo "$(COLOR_GREEN)✓ Installation complete!$(COLOR_RESET)"
	@echo "  Libraries: $(INSTALL_LIB_DIR)"
	@echo "  Headers:   $(INSTALL_INC_DIR)"

uninstall:
	@echo "$(COLOR_BOLD)Uninstalling fdpricing library...$(COLOR_RESET)"
	@rm -f $(INSTALL_LIB_DIR)/lib$(LIB_NAME).a
	@rm -f $(INSTALL_LIB_DIR)/lib$(LIB_NAME).so
	@rm -f $(INSTALL_INC_DIR)/fdpricing.h
	@ldconfig 2>/dev/null || true
	@echo "$(COLOR_GREEN)✓ Uninstall complete!$(COLOR_RESET)"

# ============================================================================
# Clean targets
# ============================================================================

clean:
	@echo "$(COLOR_YELLOW)Cleaning build artifacts...$(COLOR_RESET)"
	@rm -rf $(BUILD_DIR)
	@echo "$(COLOR_GREEN)✓ Clean complete!$(COLOR_RESET)"

distclean: clean
	@echo "$(COLOR_YELLOW)Removing all generated files...$(COLOR_RESET)"
	@rm -f *~ src/*~ include/*~ examples/*~ tests/*~
	@rm -f .*.swp src/.*.swp include/.*.swp examples/.*.swp tests/.*.swp
	@echo "$(COLOR_GREEN)✓ Distclean complete!$(COLOR_RESET)"

# ============================================================================
# Development helpers
# ============================================================================

.PHONY: print-vars format check

# Print Makefile variables (debugging)
print-vars:
	@echo "CC:            $(CC)"
	@echo "BUILD_MODE:    $(BUILD_MODE)"
	@echo "CFLAGS:        $(CFLAGS)"
	@echo "LDFLAGS:       $(LDFLAGS)"
	@echo "TEST_LDFLAGS:  $(TEST_LDFLAGS)"
	@echo "CORE_SRCS:     $(CORE_SRCS)"
	@echo "ALL_OBJS:      $(ALL_OBJS)"
	@echo "TEST_SRCS:     $(TEST_SRCS)"
	@echo "TEST_BINS:     $(TEST_BINS)"
	@echo "EXAMPLE_BINS:  $(EXAMPLE_BINS)"

# Format code (if clang-format is available)
format:
	@if command -v clang-format >/dev/null 2>&1; then \
		echo "$(COLOR_CYAN)Formatting code...$(COLOR_RESET)"; \
		find src include examples tests -name '*.c' -o -name '*.h' | xargs clang-format -i; \
		echo "$(COLOR_GREEN)✓ Format complete!$(COLOR_RESET)"; \
	else \
		echo "$(COLOR_RED)clang-format not found$(COLOR_RESET)"; \
	fi

# Static analysis (if available)
check:
	@if command -v cppcheck >/dev/null 2>&1; then \
		echo "$(COLOR_CYAN)Running static analysis...$(COLOR_RESET)"; \
		cppcheck --enable=all --suppress=missingIncludeSystem -I include src/; \
		echo "$(COLOR_GREEN)✓ Check complete!$(COLOR_RESET)"; \
	else \
		echo "$(COLOR_YELLOW)cppcheck not found, skipping static analysis$(COLOR_RESET)"; \
	fi
