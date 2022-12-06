CC :=gcc
CFLAGS := -fPIC -shared -Wall -Wextra 
SRC := src/C
OBJ := src/C
ifeq ($(OS),Windows_NT)
	EXT := .dll
else
	EXT := .so
endif
SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%$(EXT), $(SOURCES))

all: $(OBJECTS)
	
$(OBJ)/%$(EXT): $(SRC)/%.c
	$(CC) $(CFLAGS) -o $@  -c $<

clean: 
	rm -f $(OBJECTS)
