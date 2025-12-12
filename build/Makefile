# Compiler and flags
CC      = gcc
CFLAGS  = -Wall -Iinclude -g -Wextra
LDFLAGS = -lm -lgmp

# Directories
SRC_DIR = ../src
OBJ_DIR = obj
BIN_DIR = bin

# Files
SRC     = $(wildcard $(SRC_DIR)/*.c)
OBJ     = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))
TARGET  = $(BIN_DIR)/a.out

# Default target
all: $(TARGET)

# Link object files
$(TARGET): $(OBJ) | $(BIN_DIR)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create directories if they don't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Run program
run: $(TARGET)
	./$(TARGET)

# Clean build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean run

