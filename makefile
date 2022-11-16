CC :=gcc
LDFLAGS := -fPIC -shared
C_SOURCES :=$(shell find . -name '*.c')
C_EXECUTABLE :=$(C_SOURCES:.c=)
 
all:$(C_EXECUTABLE)

$(C_EXECUTABLE):$(C_SOURCES)
		$(CC) $< $(LDFLAGS) $(CFLAGS) -o $@

clean:
		rm -rf $(EXECUTABLE)