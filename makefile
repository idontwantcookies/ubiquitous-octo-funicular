CC = g++
FLAGS = -lgmp -lgmpxx -Wall -pedantic -g3

all:
	$(CC) -o generator.out src/generator.cpp $(FLAGS)

test: src/*.py
	coverage run -m pytest
	coverage html

lint: **/*.py
	pylint **/*.py
