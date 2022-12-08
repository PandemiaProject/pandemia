CC :=gcc
CFLAGS := -fPIC -shared -Wall -Wextra 
SRC := src/pandemia/C/src
OBJ := src/pandemia/C/build

ifeq ($(OS),Windows_NT)
	EXT := .dll
else
	EXT := .so
endif
SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%$(EXT), $(SOURCES))

all: $(OBJECTS)
	
$(OBJ)/%$(EXT): $(SRC)/%.c
	mkdir -p $(OBJ)

	$(CC) $(CFLAGS) -o $@  $<

clean: 
	rm -f $(OBJECTS)
