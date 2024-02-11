CC = clang
CFLAGS = -Wall -pedantic -std=c99

all: phylib

phylib.o: phylib.c phylib.h
	$(CC) $(CFLAGS) -c phylib.c -fPIC -o phylib.o

libphylib.so: phylib.o
	$(CC) $(CFLAGS) phylib.o -shared -o libphylib.so


clean: 
	rm -if *.o *.so phylib gphylib

phylib: phylib.o
	$(CC) $(CFLAGS) phylib.o -o phylib -lm
