/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 9                *
*                                                          *
************************************************************
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <fstream>
#include <math.h>

using namespace std;

#define Vector array<long double, 65>

void fillXY(Vector &x, Vector &y);
long double newtonSymbol(int x, int y);
void calculateW(Vector &w, const int d);
void calculateResult(const Vector &x, const Vector &y, const Vector &w);
long double factorial(int x);
void writeToFile(const Vector &x, const Vector &y, const char *f);

int main()
{
    int d = 3;

    Vector x;      // x-y funkcji zmieniające się o 1/32
    Vector y;      // wartości funkcji 1/(1+5*x*x)
    Vector w = {}; // wektor wag

    // funkcja interpolowana
    fillXY(x, y);
    writeToFile(x, y, "Zadanie9_function.txt");

    // wielomian interpolacyjny
    calculateW(w, d);
    calculateResult(x, y, w);

    return 0;
}

// wypełnia wektory wartościami x i wartościami funkcji
void fillXY(Vector &x, Vector &y)
{
    int i = 0;

    for (long double j = -1; j <= 1; j = j + 1.0 / 32)
    {
        x[i] = j;
        y[i] = (1 / (1 + 5 * j * j));
        i++;
    }
}

// zapisuje wektor x i wartości funkcji y do pliku
void writeToFile(const Vector &x, const Vector &y, const char *f)
{
    ofstream file;
    file.open(f);

    for (int i = 0; i < x.size(); i++)
        file << x[i] << "       " << y[i] << endl;

    file.close();
}

// liczy wagi do algorytmu
void calculateW(Vector &w, const int d)
{
    long double sum;

    for (int k = 0; k < w.size(); k++)
    {
        w[k] = pow(-1, k - d);
        sum = 0;

        for (int i = k - d; i <= k; i++)
        {
            if (i < 0)
                continue;
            if (i == w.size() - d)
                break;

            sum += newtonSymbol(d, k - i);
        }

        w[k] = w[k] * sum;
    }
}

// oblicza symbol newtona
long double newtonSymbol(int x, int y)
{
    return factorial(x) / (factorial(y) * factorial(x - y));
}

// oblicza silnie
long double factorial(int x)
{
    if (x == 0)
        return 1;

    if (x < 2)
        return x;

    return x * factorial(x - 1);
}

// oblicza wartości y i zapisuje je do pliku
void calculateResult(const Vector &x, const Vector &y, const Vector &w)
{
    long double upperSum;
    long double lowerSum;
    long double temp;

    ofstream file;
    file.open("Zadanie9_splajn.txt");

    for (long double i = -1; i <= 1; i = i + 0.01)
    {
        upperSum = 0;
        lowerSum = 0;
        temp = 0;

        for (int k = 0; k < w.size(); k++)
        {
            temp = w[k] / (i - x[k]);
            upperSum += temp * y[k];
            lowerSum += temp;
        }

        file << i << "      " << upperSum / lowerSum << endl;
    }

    file.close();
}