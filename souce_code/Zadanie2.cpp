/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 2                *
*                                                          *
************************************************************
*/

#include <iostream>
#include <vector>  // vector
#include <fstream> // fstream
#include <math.h>  // sqrt
#include <numeric> // inner_product
#include <chrono>  // chrono::high_resolution_clock
#include <stdio.h>
#include <stdlib.h>
#include <iomanip> // zawiera funkcję setprecision()

using namespace std;

#define SIZE 128
#define Vector vector<long double>
#define Matrix vector<vector<long double>>

void fill(Matrix &A, Vector &x, Vector &b);
void GaussSeidel(const Matrix &A, Vector &x, const Vector &b, const int no_of_iters);
void ConjugateGradients(const Matrix &A, Vector &x, const Vector &b, const int no_of_iters);
long double sum(int i, const Matrix &A, const Vector &x);
long double normCalculation(const Vector &prev_x, const Vector &x);
void writeToFile(const Vector &x, const char *f);
long double scalarProduct(const Vector &a, const Vector &b);

int main(int argc, char const *argv[])
{
    // kontrola liczby argumentów progamu
    if (argc != 2)
    {
        cout << "Niepoprawna liczba argumentów wywołania programu!" << endl;
        cout << "Prawidlowe wywolanie programu: ./Zadanie2 <liczba iteracji>" << endl;
        exit(-1);
    }

    Matrix A(SIZE);

    for (int i = 0; i < SIZE; i++)
        A[i].resize(SIZE); // zmiana wielkości

    Vector x(SIZE);
    Vector b(SIZE);

    cout << "Czas wykonania dla metod: " << endl;

    // wypełnianie wektorów
    fill(A, x, b);

    auto start = chrono::high_resolution_clock::now();
    GaussSeidel(A, x, b, stoi(argv[1]));
    auto end = chrono::high_resolution_clock::now();

    cout << "a) Gaussa-Seidela:"            << "          " << chrono::duration<double>(end-start).count() << endl;

    // wypełnianie od nowa dla drugiej metody
    fill(A, x, b);

    start = chrono::high_resolution_clock::now();
    ConjugateGradients(A, x, b, stoi(argv[1]));
    end = chrono::high_resolution_clock::now();

    cout << "b) Gradientów Sprzezonych:"    << "  " << chrono::duration<double>(end-start).count() << endl;

    return 0;
}

void fill(Matrix &A, Vector &x, Vector &b)
{
    for (int i = 0; i < x.size(); i++)
    {
        x[i] = -1;
        b[i] = 1;

        if (i == 0)
            A[i][i + 1] = 1;

        else if (i == x.size() - 1)
            A[i][i - 1] = 1;

        else
        {
            A[i][i - 1] = 1;
            A[i][i + 1] = 1;
        }

        A[i][i] = 4;

        if (i < x.size() - 4)
        {
            A[i][i + 4] = 1;
            A[i + 4][i] = 1;
        }
    }
}

// metoda Gaussa-Seidela dla <no_of_iters> iteracji
void GaussSeidel(const Matrix &A, Vector &x, const Vector &b, const int no_of_iters)
{
    int iter_num = 0; // ile już było iteracji
    Vector prev_x;    // wektor poprzednich wartości x
    Vector norms;     // wektor dla obliczonych norm
    prev_x.resize(x.size());

    for (int i = 0; i < x.size(); i++) // dodawanie nowego elementu na koniec wektora
        prev_x.push_back(x[i]);

    for (int it = 0; it < no_of_iters; it++)
    {
        // obliczamy kolejną iterację x
        for (int i = 0; i < x.size(); i++)
        {
            prev_x[i] = x[i];
            x[i] = (b[i] - sum(i, A, x)) / A[i][i];
        }

        // każdorazowo obliczamy normę
        norms.push_back(normCalculation(prev_x, x));
    }

    // zapisujemy wektor norm do pliku (do późniejszego wykresu)
    writeToFile(norms, "Zadanie2_GaussSeidel.txt");
}

// suma składników mnożenia macierzy i wektora
long double sum(int i, const Matrix &A, const Vector &x)
{
    long double sum = 0;

    if (i == 0)
        sum = +A[0][1] * x[1] + A[0][4] * x[4];

    else if (i < 4)
        sum = A[i][i - 1] * x[i - 1] + A[i][i + 1] * x[i + 1] + A[i][i + 4] * x[i + 4];

    else if (i < A.size() - 4)
        sum = A[i][i - 4] * x[i - 4] + A[i][i - 1] * x[i - 1] + A[i][i + 1] * x[i + 1] + A[i][i + 4] * x[i + 4];

    else if (i < A.size() - 1)
        sum = A[i][i - 4] * x[i - 4] + A[i][i - 1] * x[i - 1] + A[i][i + 1] * x[i + 1];

    else
        sum = A[i][i - 4] * x[i - 4] + A[i][i - 1] * x[i - 1];

    return sum;
}

//  norma różnicy wektorów
long double normCalculation(const Vector &prev_x, const Vector &x)
{
    long double norm;
    long double res = 0;

    for (int i = 0; i < x.size(); i++)
    {
        norm = x[i] - prev_x[i];
        res += norm * norm;
    }

    res = sqrt(res);
    return res;
}

// zapis do pliku (do stworzenia wykresu)
void writeToFile(const Vector &x, const char *f)
{
    ofstream file;

    file.open(f);

    for (int i = 0; i < x.size(); i++)
        file << i + 1 << ". " << x[i] << endl;

    file.close();
}

// metoda Gradientów Sprzężonych wykonana <no_of_iters> razy
void ConjugateGradients(const Matrix &A, Vector &x, const Vector &b, const int no_of_iters)
{
    Vector r = b;
    Vector p = r;
    Vector prev_x = x; //wektor poprzednich wartosci x
    Vector norms;      //wektor dla obliczonych norm
    Vector Apk(x.size());

    for (int i = 0; i < no_of_iters; i++)
    {
        Vector prev_r = r;

        // obliczanie Apk
        for (int j = 0; j < x.size(); j++)
            Apk[j] = scalarProduct(A[j], p);

        long double alfa = scalarProduct(r, r) / scalarProduct(p, Apk);

        // następny wynik x
        for (int j = 0; j < x.size(); j++)
        {
            prev_x[j] = x[j];
            x[j] = x[j] + alfa * p[j];
        }

        for (int j = 0; j < x.size(); j++)
        {
            r[j] = r[j] - alfa * Apk[j];
        }

        long double beta = scalarProduct(r, r) / scalarProduct(prev_r, prev_r);

        //kolejne p
        for (int j = 0; j < x.size(); j++)
        {
            p[j] = r[j] + beta * p[j];
        }

        //obliczenie normy
        norms.push_back(normCalculation(prev_x, x));
    }

    // zapisujemy wektor norm do pliku (do późniejszego wykresu)
    writeToFile(norms, "Zadanie2_ConjGrad.txt");
}

// oblicza iloczyn skalarny wektorow a i b
long double scalarProduct(const Vector &a, const Vector &b)
{
    return inner_product(a.begin(), a.end(), b.begin(), (long double)0.0);
}