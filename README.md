# Problema do logaritmo discreto

## Desenvolvimento do trabalho

<!-- TODO -->

## Módulos

O programa é constituído dos seguintes módulos:

- `base`: Possui funções aritméticas simples, para o cálculo de raiz quadrada, logaritmo na base 10 e MDC, por exemplo.
- `util`: funções úteis para medir tempo e para sair do programa em caso de erro / timeout.
- `mod`: aqui situam-se funções fundamentais relativas à aritmética modular, como inverso modular, potenciação modular, subgrupo, etc. Além disso, aqui se encontra a função `find_generator` usada pelo programa na segunda etapa do algoritmo. O algoritmo estendido de Euclides (do módulo `base`) é usado com frequência nesse módulo.
- `primality`: funções relacionadas a testes de primalidade, incluindo o teste de Miller-Rabin e o crivo de Eratóstenes. O teste de Miller-Rabin usa o módulo `base` para calcular log na base 10, e o crivo de Eratóstenes usa o mesmo módulo para calcular a raiz quadrada inteira.
- `factorization`: como o próprio nome já diz, aqui se encontram as funções responsáveis por fatorar um número inteiro, incluindo a função totiente. A fatoração de n - 1, onde n é o primo encontrado na primeira etapa, é utilizada em diversos algoritmos subsequentes para acelerar os cálculos. A fatorização usa o teste de primalidade Miller-Rabin para checar se o número foi completamente fatorado.
- `discrete_log`: algoritmos baby-step, giant-step e Pohlig-Hellman, usados para resolver o problema do logaritmo discreto na terceira etapa do programa. O cálculo do logaritmo discreto depende da fatoração feita pelo módulo `factorization`, além das funções de potenciação modular e a solução de congruências pelo Teorema Chinês do Resto, no módulo `mod`.
- `main`: o script principal que recebe como entrada (N, h) para encontrar:
    - o menor n primo tal que n > N;
    - a fatoração de n - 1;
    - um gerador g de Zn, ou um elemento com ordem alta;
    - o logaritmo discreto de h na base g mod n.

## Formato de entrada e como utilizar o programa
A entrada deve ser dada por dois números inteiros, N e a, separados por _whitespace_. Pode-se rodar o programa, e depois fornecer os inteiros sequencialmente, ou passar um arquivo contendo os 2 inteiros separados por _whitespace_, usando `cat [nome do arquivo].txt | python main.py`, por exemplo. A pasta `test/` possui 3 exemplos de entrada formatados em arquivos válidos: `in1.txt`, `in2.txt` e `in3.txt`.

Também é possível invocar o programa diretamente e inserir os valores manualmente, através de da linha de comando, com o comando `python main.py`.

```
Insira o valor de N: 90
Insira o valor de a: 4
```

## Formato de saída
A saída é dividida em 4 partes:

1. O cálculo do menor primo n maior que N, seguido do número de repetições feitas do teste de Miller, e, por fim, o tempo de execução para encontrar n.

```
Menor primo n maior que N: 97
Repetições de Miller-Rabin usadas: 10
Tempo de execução: 0.057ms.
```

2. A decomposição de n - 1 em primos, seguida do tempo gasto para seu cálculo.

```
Decomposição de n - 1 em primos: Counter({2: 5, 3: 1})
Tempo de execução: 0.124ms.
```

3. Um gerador g para o grupo Zn, seguido do tempo gasto para encontrá-lo.

```
Gerador: g= 58
Tempo de execução: 0.015ms.
```

Além disso, caso não seja possível encontrar um gerador em menos de 15 segundos, a saída será um inteiro g tal que g tem a maior ordem que o programa for capaz de calcular. O programa encerra o cálculo aqui caso não consiga encontrar um gerador, pois, sem gerador, nem sempre é possível encontrar o logaritmo discreto de h.

4. O log discreto de h na base g mod n, seguido do tempo de exeução gasto. Novamente, aqui há um tempo limite de 15 segundos, e, caso não seja possível calcular o log, o programa se encerra.

```
Log discreto de h na base g: 28
Tempo de execução: 0.099ms.
```

## Complexidade

As 4 funções principais executadas são: 
- `prime_miller_rabin()`, que tem complexidade O(rep * log²(n)). Como adotamos rep = log(n), a complexidade total é O(log³(n));
- `factors()`, de complexidade O(log(n) * sqrt(p)), onde p é o maior fator primo de n - 1;
- `find_generator()`, de complexidade O(#########) asdfdas
- `pohlig_hellman()`, de complexidade O(r * sqrt(p)).

<!-- TODO: Complexidade de find_generator() -->


## Listagem do programa-fonte

<!-- TODO -->