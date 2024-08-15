# Problema do logaritmo discreto - TP1 de Álgebra A

## Autores

- Chrystian Paulo Ferreira de Melo
- Felipe Ribas Muniz
- João Marcos Rezende
- Sanny Cristiane Moreira de Sales

## Desenvolvimento do trabalho

O trabalho foi desenvolvido em python 3.10. Para gerar os casos de teste do código, usamos a biblioteca GMP de C++ para determinar primos grandes, potências mod n, inversos modulares e testes de primalidade. O arquivo generator.cpp contém o código fonte responsável pelos casos de teste gerados, que foram então usados em python para testar o código do aplicativo de fato.

Diversos casos de teste foram criados no arquivo `tp_test.py` em casos onde sabíamos a solução, geralmente para números menores do que aqueles gerados pelo GMP e em que foi possível fazer os cálculos na mão. A base do trabalho foi feita sobre o algoritmo estendido de Euclides para calcular o MDC e sobre o algoritmo de exponenciação binária. Com eles, somos capazes de gerar grupos em Zp, e de calcular inversos modulares. Além disso, o teste de primalidade de Miller-Rabin foi usado para encontrar um p grande.

Uma vez que essas 3 bases foram construídas, temos um p e um grupo Zp, onde todo elemento tem inverso multiplicativo. A partir dessa premissa, trabalhamos com o totiente de p, a fatorização do totiente em potências de primos usando o algoritmo Pohlard's rho, buscamos um gerador g para Zp, e, por fim, usamos Baby-Step, Giant-Step, combinado com Pohlig-Hellman, para calcular o logaritmo discreto de qualquer elemento h ∈ Zp na base g.

Os testes foram usados com frequência para refatorar o código e realizar mudanças necessárias sempre que algum erro de lógica era encontrado, ou quando era necessário aprimorar a eficiência de um algoritmo. Eles podem não ser exaustivos; além disso, como alguns desses algoritmos são probabilísticos, os testes são instáveis, embora eles dificilmente falhem duas vezes consecutivas.

## Módulos

O programa é constituído dos seguintes módulos:

### base

Possui funções aritméticas simples, para o cálculo de raiz quadrada, logaritmo na base 10 e MDC, por exemplo.

### util

Funções úteis para medir tempo e para sair do programa em caso de erro / timeout.

### modular_arithmetic

Aqui situam-se funções fundamentais relativas à aritmética modular, como inverso modular, potenciação modular, subgrupo, etc. Além disso, aqui se encontra a função find_generator usada pelo programa na segunda etapa do algoritmo. O algoritmo estendido de Euclides (do módulo base) é usado com frequência nesse módulo.

### primality

Funções relacionadas a testes de primalidade, incluindo o teste de Miller-Rabin e o crivo de Eratóstenes. O teste de Miller-Rabin usa o módulo base para calcular log na base 10, e o crivo de Eratóstenes usa o mesmo módulo para calcular a raiz quadrada inteira.

### factorization

Como o próprio nome já diz, aqui se encontram as funções responsáveis por fatorar um número inteiro, incluindo a função totiente. A fatoração de n - 1, onde n é o primo encontrado na primeira etapa, é utilizada em diversos algoritmos subsequentes para acelerar os cálculos. A fatorização usa o teste de primalidade Miller-Rabin para checar se o número foi completamente fatorado.

### discrete_log 

Algoritmos baby-step, giant-step e Pohlig-Hellman, usados para resolver o problema do logaritmo discreto na terceira etapa do programa. O cálculo do logaritmo discreto depende da fatoração feita pelo módulo factorization, além das funções de potenciação modular e a solução de congruências pelo Teorema Chinês do Resto, no módulo modular_arithmetic.

### linalg

Operações em vetores e matrizes, como multiplicação de matrizes, redução à forma escalonada (RREF) e cálculo do kernel de uma matriz.

### quadratic_sieve

Implementação do crivo quadrático para fatoração de números grandes, usando resíduos quadráticos e álgebra linear para encontrar fatores não triviais.

### rsa

Sistema de criptografia RSA, incluindo funções para gerar chaves públicas e privadas e codificar e decodificar mensagens usando aritmética modular.
### tp1

Script principal que recebe como entrada (N, h) para encontrar:
    - o menor n primo tal que n > N;
    - a fatoração de n - 1;
    - um gerador g de Zn, ou um elemento com ordem alta;
    - o logaritmo discreto de h na base g mod n.

### tp2
Script principal para fatorar um número inteiro usando o crivo quadrático.

Recebe como entrada o número n, e retorna um fator desse número.

## Como utilizar o programa

É necessário ter o python 3.10 instalado localmente. Para rodar os testes, também é preciso ter a biblioteca pytest instalada.


Para executar o programa, entre na pasta-raiz do projeto e rode o comando

```bash
cd {nome-do-diretorio}
python tp1.py
```

Para checar se tudo está nos conformes, basta trocar para o diretório-raiz do projeto e executar o pytest.

```bash
cd {nome-do-diretorio}
pytest
```
Para alimentar as informações de entrada para o aplicativo, é possível usar a pipe do bash para passar os argumentos para o programa.

```bash
cat test/in3.txt | python tp1.py
```

## Formato de entrada

A entrada deve ser dada por dois números inteiros, N e a, separados por _whitespace_. Pode-se rodar o programa, e depois fornecer os inteiros sequencialmente, ou passar um arquivo contendo os 2 inteiros separados por _whitespace_, usando `cat [nome do arquivo].txt | python tp1.py`, por exemplo. A pasta `test/` possui 3 exemplos de entrada formatados em arquivos válidos: `in1.txt`, `in2.txt` e `in3.txt`.

Também é possível invocar o programa diretamente e inserir os valores manualmente, através de da linha de comando, com o comando `python tp1.py`.

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
- `find_generator()`, de complexidade O(log²(p))
- `pohlig_hellman()`, de complexidade O(r * sqrt(p)), onde r é o número de fatores primos de n - 1 e p é o maior fator primo de n - 1.


## Listagem do programa-fonte

Os arquivos do código-fonte são:

- src/ base.py, discrete_log.py, factorization.py, linalg.py, modular_arithmetic.py, primality.py, quadratic_sieve.py, rsa.py, util.py: módulos essenciais da aplicação.
- tp2.py: código principal da aplicação.
- tests/*.py: casos de teste dos módulos essenciais.
- generator.cpp: código-fonte para gerar primos grandes para os casos de teste
- big_numbers/*.txt: listas de testes usadas por tp_test.py e para os testes manuais.
