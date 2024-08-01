CC = g++
FLAGS = -lgmp -lgmpxx -Wall -pedantic -g3

all:
	$(CC) -o generator.out generator.cpp $(FLAGS)

test: *.py
	coverage run -m pytest
	coverage html

lint: *.py
	pylint *.py
