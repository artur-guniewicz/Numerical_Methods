/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 11               *
*                                                          *
************************************************************
*/

#include <iostream>
#include <math.h>  // zawiera funkcje pow, sqrt, sin, abs
#include <vector>  // umożliwia tworzenie wektorów
#include <iomanip> // zawiera funkcję setprecision
#include <stdio.h>
#include <stdlib.h>

#define e 2.7182818284590452354
#define pi 3.14159265358979323846

using namespace std;

long double trapezoidArea(long double, long double, int);
long double function(long double);

int main(int argc, char const *argv[])
{
    // kontrola liczby argumentów progamu
    if (argc != 4)
    {
        cout << "Niepoprawna liczba argumentów wywołania programu!" << endl;
        cout << "Prawidlowe wywolanie programu: ./Zadanie11 <liczba iteracji> <dolna granica całki> <przybliżenie 10^x>" << endl;
        exit(-1);
    }

    // dokładność 10^<argument wywolania>
    int precision = stoi(argv[3]);
    long double tolerance = pow(10, precision);

    // szukanie A takiego, że e^-A < 10^-7
    int i = -1;
    long double rest = 0;
    while (true)
    {
        rest = pow(e, i);

        if (rest < tolerance)
            break;

        i--;
    }

    cout << "A = " << -i << endl;

    // górna granica całki
    long double upper = -i;

    // dolna granica całki
    long double lower = stoi(argv[2]);

    // limit iteracji
    long double limitOfIters = stoi(argv[1]);

    vector<long double> rombergs;      // obliczanie przybliżeń do obecnego wiersza
    vector<long double> prev_rombergs; // zapamiętuje poprzednie przybliżenia (poprzedni wiersz)

    // ustawianie precyzji
    cout << setprecision(-precision) << fixed;

    // oblicz pierwszy wyraz
    rombergs.push_back(trapezoidArea(lower, upper, 0));

    cout << "k = 0 -> " << rombergs[0] << endl;

    int tempK;
    long double temp;
    int k;

    for (k = 1; k < limitOfIters - 1; k++)
    {
        // zapamiętanie poprzedniego przybliżenia
        prev_rombergs = rombergs;

        // czyścimy obecny wektor
        rombergs.clear();
        cout << "k = " << k << " -> ";
        tempK = k;

        for (int n = 0; n <= k; n++)
        {
            if (n == 0)
            {
                rombergs.push_back(trapezoidArea(lower, upper, k));
                cout << rombergs[0] << " ";
                tempK--;
            }
            else
            {
                temp = (pow(4, n) * rombergs[n - 1] - prev_rombergs[n - 1]) / (pow(4, n) - 1);
                cout << temp << " ";
                rombergs.push_back(temp);
                tempK--;
            }
        }

        cout << endl;

        if (abs(rombergs[k] - prev_rombergs[k - 1]) < tolerance)
        {
            cout << endl << "Pole całki w przedziale [" << lower << "," << upper << "] wynosi: " << rombergs[k] << " (z dokładnością do " << tolerance << ")" << endl;
            break;
        }
    }

    long double result = rombergs[k] + rest;
    cout << endl << "Pole łaczne całki w przedziale [0,inf) wynosi: " << result << " (z dokładnością do " << tolerance << ")" << endl;
    return 0;
}

// korzysta ze złożonej metody trapezow by podzielic całkę na k przedziałów
long double trapezoidArea(long double low, long double high, int k)
{
    // ilość punktów podziału
    long double N = pow(2, k);
    long double area = 0;

    // przedział pomiędzy argumentami funkcji
    long double h = (high - low) / N;

    for (long double i = low; i <= high; i = i + h)
    {
        if (i == low || i == high)
            area += function(i) / 2;

        else
            area += function(i);
    }

    return area * h;
}

// oblicza wartość funkcji podcałkowej dla danego x
long double function(long double x)
{
    return sin(pi * (1 + sqrt(x)) / (1 + pow(x, 2))) * pow(e, -x);
}