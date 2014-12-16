SHELL = /bin/bash
CC    = gcc
CFLAGS  = -Wall

OPTRUNDBG = -d 1 -l 15 -k 2 -c input/sequences_test_1.fasta -o output/resultats_test.txt -i 20

TEST_1_1 = -d 2 -l 6 -k 3 -c input/sequences_test_1.fasta -o output/resultats_test_1_1.txt -i 10

TEST_1_2 = -d 2 -l 6 -k 1 -c input/sequences_test_1.fasta -o output/resultats_test_1_2.txt -i 10

TEST_2 = -d 1 -l 4 -k 3 -c input/sequences_test_2.fasta -o output/resultats_test_2.txt -i 1000

TEST_3 = -d 2 -l 7 -k 3 -c input/sequences_test_3.fasta -o output/resultats_test_3.txt -i 10

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
	$(CC)  $(CFLAGS) $(TARGETC) $(OBJECTS) -o $(TARGET) -lm -g

$(OBJDIR)/%.o: $(FCTDIR)/%.c $(HEADERS) 
	$(CC)  $(CFLAGS) -c -o $@ $< -g

clean:
	rm $(OBJDIR)/*.o -rf

mrproper: clean
	rm -rf $(TARGET)

dbg: $(TARGET)
	valgrind --leak-check=full --track-origins=yes  --show-leak-kinds=all ./$(TARGET) $(OPTRUNDBG)

time: $(TARGET)
	time ./$(TARGET) $(OPTRUNDBG)

test_1_1: $(TARGET)
	valgrind --leak-check=full --track-origins=yes  --show-leak-kinds=all ./$(TARGET) $(TEST_1_1)

test_1_2: $(TARGET)
	valgrind --leak-check=full --track-origins=yes  --show-leak-kinds=all ./$(TARGET) $(TEST_1_2)

test_2: $(TARGET)
	valgrind --leak-check=full --track-origins=yes  --show-leak-kinds=all ./$(TARGET) $(TEST_2)

test_3: $(TARGET)
	valgrind --leak-check=full --track-origins=yes  --show-leak-kinds=all ./$(TARGET) $(TEST_3)

.PHONY: clean rundbg mrproper
