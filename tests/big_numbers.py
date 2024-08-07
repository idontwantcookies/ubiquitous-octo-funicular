extgcd = []
with open("test/extended_gcd.txt") as file:
    for line in file:
        extgcd.append([int(x) for x in line.split()])

bsgs = []
with open("test/bsgs.txt") as file:
    for line in file:
        bsgs.append([int(x) for x in line.split()])

modular_inverse = []
with open("test/inverso_modular.txt") as file:
    for line in file:
        modular_inverse.append([int(x) for x in line.split()])

fastexp = []
with open("test/exp_binaria.txt") as file:
    for line in file:
        fastexp.append([int(x) for x in line.split()])

primes = []
with open("test/primes.txt") as file:
    for line in file:
        primes.append([int(x) for x in line.split()])
