CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

# Specify the target executable and the source files needed to build it
my_app: symnmf.o symnmf.h
	$(CC) -o my_app symnmf.o $(CFLAGS)

# Specify the object files that are generated from the corresponding source files
symnmf.o: symnmf.c
	$(CC) -c symnmf.c $(CFLAGS)
