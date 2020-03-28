/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 1a               *
*                                                          *
************************************************************
*/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip> // zawiera funkcję setprecision()

using namespace std;

// drukuje równanie Ax=b
void printEquation(const long double A[7][7], const long double x[7], const long double b[7])
{
    for (int i = 0; i < 7; i++)
    {
        cout << setprecision(3) << fixed;

        cout << "[ ";

        for (int j = 0; j < 7; j++)
            cout << A[i][j] << " ";

        cout << "] [ " << x[i] << " ] [ " << b[i] << " ]" << endl;
    }

    cout << endl;
}

// drukuje wektor
void printVector(string name, const long double x[7])
{
    cout << setprecision(10);

    cout << "Wektor " << name << ":" << endl;

    cout << name << " = [ ";
    for (int i = 0; i < 7; i++)
        cout << x[i] << " ";
    cout << "]" << endl;
}

// faktoryzacja Choleskiego
void CholeskyDecomposition(const long double A[7][7], long double *x, const long double *b)
{
    long double C[7][7] = {0};

    C[0][0] = sqrt((A[0][0]));

    for (int i = 1; i < 7; i++)
    {
        C[i][i - 1] = A[i][i - 1] / C[i - 1][i - 1];
        C[i][i] = sqrt((A[i][i] - C[i][i - 1] * C[i][i - 1]));
    }

    double long y[7] = {0};

    //forwardsubstitution Cy=b, gdzie y=C^Tx
    y[0] = b[0] / C[0][0];

    for (int i = 1; i < 7; i++)
        y[i] = (b[i] - y[i - 1] * C[i][i - 1]) / C[i][i];

    //backsubstitution C^Tx=y
    x[6] = y[6] / C[6][6];

    for (int i = 5; i >= 0; i--)
        x[i] = (y[i] - C[i + 1][i] * x[i + 1]) / C[i][i];
}

int main()
{
    long double A[7][7] = {0};
    long double x[7] = {0};
    long double b[7] = {0};

    // uzupełnienie macierzy A i b
    for (int i = 0; i < 7; i++)
    {
        if (i == 0)
            A[i][i + 1] = 1;

        else if (i == 6)
            A[i][i - 1] = 1;

        else
        {
            A[i][i - 1] = 1;
            A[i][i + 1] = 1;
        }

        A[i][i] = 4;
        b[i] = i + 1;
    }

    cout << endl;

    CholeskyDecomposition(A, x, b);

    printEquation(A, x, b);

    printVector("x", x);

    return 0;
}