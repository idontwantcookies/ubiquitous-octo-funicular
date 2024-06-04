CC = g++
FLAGS = -lgmp -lgmpxx -Wall -pedantic -g3
all:
	$(CC) -o generator.out generator.cpp $(FLAGS)
