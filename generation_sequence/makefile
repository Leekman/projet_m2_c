CC=gcc
CFLAGS=-Wall 
LDFLAGS=
EXEC=generateur_seq
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)




all: $(EXEC)

main.o: main.c 

generateur_seq: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ 

%.o: %.c
	$(CC) $(CFLAGS) -o opt/$@ -c fct/$< 


.PHONY: clean mrproper

clean:
	rm -rf opt/*.o

mrproper: clean
	rm -rf bin/generateur_seq