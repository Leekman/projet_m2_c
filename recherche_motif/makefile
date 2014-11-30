SHELL = /bin/bash
CC    = gcc
CFLAGS  = -Wall

OPTRUNDBG = -l 4 -k 2 -c input/sequences.fasta

BINDIR = bin
LIBDIR = lib
SRCDIR = src
OBJDIR = obj
FCTDIR = fct
NAME = $(shell basename $(CURDIR))

TARGET  = $(BINDIR)/$(NAME)
TARGETC = $(SRCDIR)/$(NAME).c 
SOURCES=$(wildcard $(FCTDIR)/**/*.c $(FCTDIR)/*.c)
OBJECTS=$(patsubst $(FCTDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))
HEADERS = $(wildcard $(LIBDIR)/**/*.h $(LIBDIR)/*.h)

all: $(TARGET)

$(TARGET): $(TARGETC) $(OBJECTS) 
	$(CC)  $(CFLAGS) $(TARGETC) $(OBJECTS) -o $(TARGET) -lm

$(OBJDIR)/%.o: $(FCTDIR)/%.c $(HEADERS) 
	$(CC)  $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJDIR)/*.o -rf

mrproper: clean
	rm -rf $(TARGET)

dbg: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(OPTRUNDBG)

time: $(TARGET)
	time ./$(TARGET) $(OPTRUNDBG)

.PHONY: clean rundbg mrproper
