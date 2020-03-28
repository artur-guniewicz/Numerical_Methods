/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 8                *
*                                                          *
************************************************************
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <fstream>

using namespace std;

#define Vector array<long double, 65>

long double calcFunction(long double x);
long double base(const Vector &x, int i, long double y);
void calcLagrange(const Vector &x, const Vector &y, Vector &a);
long double calculate(const Vector &x, const Vector &a, long double y);
void writeToFile(const Vector &x, const Vector &y, const char *f);

int main()
{
    Vector x;
    Vector y;
    Vector a;

    // wypełnianie wektora punktów
    int i = 0;
    for (long double j = -1; j <= 1; j = j + (long double)1.0 / 32)
    {
        x[i] = j;
        y[i] = calcFunction(j);
        i++;
    }

    // obliczenie współczynników Lagrange'a
    calcLagrange(x, y, a);

    // zapisywanie punktow obliczonych z wielomianu interpolacyjnego
    ofstream file;
    file.open("Zadanie8_poly.txt");
    for (long double i = -1; i <= 1; i = i + 0.01)
        file << i << "      " << calculate(x, a, i) << endl;

    file.close();

    // zapiywanie wynikow do pliku
    writeToFile(x, y, "Zadanie8_function.txt");

    return 0;
}

// liczy wartosc funkcji dla danego x
long double calcFunction(long double x)
{
    long double temp = 5 * x * x + 1.0;

    return (1.0 / temp);
}

// zwraca element bazy wektora x(y), gdzie i to i-ty element bazy
long double base(const Vector &x, int i, long double y)
{
    long double base = 1;

    for (int j = 0; j < x.size(); j++)
    {
        if (j == i)
            continue;
        else
            base *= (y - x[j]);
    }

    return base;
}

// obliczenie współczynników Lagrange'a
void calcLagrange(const Vector &x, const Vector &y, Vector &a)
{
    for (int i = 0; i < a.size(); i++)
        a[i] = y[i] / base(x, i, x[i]);
}

// obliczenie wartości funkcji dla x
long double calculate(const Vector &x, const Vector &a, long double y)
{
    long double result = 0;

    for (int i = 0; i < a.size(); i++)
        result += a[i] * base(x, i, y);

    return result;
}

// zapisanie wektora x i wartości funkcji y do pliku
void writeToFile(const Vector &x, const Vector &y, const char *f)
{
    ofstream file;
    file.open(f);

    for (int i = 0; i < x.size(); i++)
        file << x[i] << "       " << y[i] << endl;
}