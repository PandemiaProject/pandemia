CC :=gcc
CFLAGS := -fPIC -shared -Wall -Wextra 
SRC := src/pandemia/cfiles
OBJ := src/pandemia/cfiles
ifeq ($(OS),Windows_NT)
	EXT := .dll
else
	EXT := .so
endif
SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%$(EXT), $(SOURCES))

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $<

$(OBJ)/%$(EXT): $(SRC)/%.c
	$(CC) -I $(SRC) -c $< -o $@

clean: 
	rm -f $(OBJECTS)
