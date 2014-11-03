SHELL = /bin/bash
CC    = gcc

FLAGS        = -std=gnu99
CFLAGS       = -Wall
DEBUGFLAGS   = -O0 -g
RELEASEFLAGS = -O2 -D NDEBUG

OPTRUNDBG = --length 9 -d 3 -k 6 -f tests/random1.fa

BINDIR = bin
LIBDIR = lib
SRCDIR = src
OBJDIR = obj
NAME = $(shell basename $(CURDIR))

TARGET  = $(BINDIR)/$(NAME)
TARGETC = $(SRCDIR)/$(NAME).c
SOURCES=$(wildcard $(LIBDIR)/**/*.c $(LIBDIR)/*.c)
OBJECTS=$(patsubst $(LIBDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))
COMMON  =
HEADERS = $(wildcard $(LIBDIR)/**/*.h $(LIBDIR)/*.h)

all: $(TARGET)

$(TARGET): $(TARGETC) $(OBJECTS)
	$(CC) $(FLAGS) $(CFLAGS) $(DEBUGFLAGS) $(TARGETC) $(OBJECTS) -o $(TARGET)

release: $(SOURCES) $(HEADERS) $(COMMON)
	$(CC) $(FLAGS) $(CFLAGS) $(RELEASEFLAGS) -o $(TARGET) $(SOURCES)

$(OBJDIR)/%.o: $(LIBDIR)/%.c $(HEADERS) $(COMMON)
	$(CC) $(FLAGS) $(CFLAGS) $(DEBUGFLAGS) -c -o $@ $<

clean:
	rm $(OBJDIR)/*.o -f

mrproper: clean
	rm $(TARGET)

rundbg: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all  ./$(TARGET) $(OPTRUNDBG)

.PHONY: clean rundbg
