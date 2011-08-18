##
##  Makefile for range-energy programs.  
##  
##  To compile:
##
##     1.  Change the value of the DEDX constant in crange.h
##     2.  make complex
##     3.  make crange
##
##  To remove stuff:
##
##     1.  make clean
## 
#CCFLAGS = -c -g -ansi -pedantic -W -Wall
CCFLAGS = -c -O3 -ansi -pedantic -W -Wall

.SUFFIXES: .c .o .r

.c.o:
	gcc $(CCFLAGS) $?


complex: complex.o
	ar ru libcomplex.a complex.o
	ranlib libcomplex.a

crange: crange.o
	gcc -o $@ crange.o libcomplex.a -lm

clean:
	rm complex.o
	rm crange.o
	rm libcomplex.a

