#include <iostream>
#include <string.h>

#include <gmpxx.h>
#include <gmp.h>

#define N 20
#define SEED 0
#define BITS 100

gmp_randclass r1(gmp_randinit_default);

using namespace std;

void gerar_mdcs()
{
    cout << "Gerando MDCs estendidos..." << endl;
    mpz_class a, b, x, y, mdc;

    for (int i = 0; i < 20; i++) {
        a = r1.get_z_bits(BITS);
        b = r1.get_z_bits(BITS);
        mpz_gcdext(mdc.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        cout << a << " " << b << " " << x << " " << y << " " << mdc << endl;
    }
}

void gerar_inverso_modular()
{
    cout << "Gerando inversos modulares..." << endl;
    mpz_class a, n, inverso, invertivel;

    for(int i=0; i<N; i++)
    {
        a = r1.get_z_bits(BITS);
        n = r1.get_z_bits(BITS);
        if(n == 0) continue;
        invertivel = mpz_invert(inverso.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
        if (!invertivel) inverso = 0;
        cout << a << " " << n << " " << inverso << endl;
    }   
}

void gerar_exponenciacao_binaria()
{
    cout << "Gerando exponenciacao binaria..." << endl;
    mpz_class b, e, n, result;

    for(int i=0; i<N; i++)
    {
        b = r1.get_z_bits(BITS);
        e = r1.get_z_bits(BITS);
        n = r1.get_z_bits(BITS);

        mpz_powm(result.get_mpz_t(), b.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
        cout << b << " " << e << " " << n << " " << result << endl;
    }
}

void gerar_primos() 
{
    clog << "Gerando números primos e compostos...\n";

    mpz_class n;
    bool primo;

    for(int i=0; i<N; i++) {
        n = r1.get_z_bits(BITS) | 1;
        primo = bool(mpz_probab_prime_p(mpz_class(n).get_mpz_t(), 20));
        cout << n << " " << primo << endl;
    }
}

int main()
{
    r1.seed(SEED);
    clog << "Gerador iniciado com seed = " <<  SEED << endl;
    gerar_mdcs();
    gerar_exponenciacao_binaria();
    gerar_inverso_modular();
    gerar_primos();
    return 0;
}
