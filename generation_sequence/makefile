SHELL = /bin/bash
CC    = gcc
CFLAGS  = -Wall

OPTRUNDBG = -d 2 -n 50 -t 100 -m TATATGCA -v 20 -o ../recherche_motif/input/sequences_generer.fasta -i output/info.txt 

TEST_1 = -d 3 -n 50 -t 50 -m ATCGTG -v 10 -o ../recherche_motif/input/sequences_test_1.fasta -i output/info_test_1.txt

TEST_2 = -d 1 -n 50 -t 50 -m AGCT -v 10 -o ../recherche_motif/input/sequences_test_2.fasta -i output/info_test_2.txt

TEST_3 = -d 2 -n 50 -t 50 -m AAAAAAA -v 10 -o ../recherche_motif/input/sequences_test_3.fasta -i output/info_test_3.txt

TEST_4 = -d 0 -n 50 -t 50 -m AAAAAAAA -v 0 -o ../recherche_motif/input/sequences_test_4.fasta -i output/info_test_4.txt

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
	$(CC)  $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJDIR)/*.o -rf

mrproper: clean
	rm -rf $(TARGET)

rundbg: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(OPTRUNDBG)

runtime: $(TARGET)
	time ./$(TARGET) $(OPTRUNDBG)

test_1: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(TEST_1)

test_2: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(TEST_2)

test_3: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(TEST_3)

test_4: $(TARGET)
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./$(TARGET) $(TEST_4)

.PHONY: clean rundbg mrproper
