/*
************************************************************
*                                                          *
*               Artur Guniewicz - Zadanie 7                *
*                                                          *
************************************************************
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>

using namespace std;

const int N = 8;

int main()
{
    int s;

    double data_x[] = {0.062500, 0.187500, 0.312500, 0.437500, 0.562500, 0.687500, 0.812500, 0.935700};
    double data_fx[] = {0.687959, 0.073443, -0.517558, -1.077264, -1.600455, -2.080815, -2.507266, -2.860307};

    gsl_vector_view x = gsl_vector_view_array(data_x, N);
    gsl_vector_view b = gsl_vector_view_array(data_fx, N);

    gsl_vector *a = gsl_vector_alloc(N);
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_matrix *V = gsl_matrix_alloc(N, N);


    // w macierzy V na pozycji (i,j) ustawiona zostaje wartość wzięta z wektora x z pozycji i
    // podniesiona do potęgi N-1-j
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            gsl_matrix_set(V, i, j, gsl_sf_pow_int(gsl_vector_get(&x.vector, i), N - 1 - j));

    gsl_linalg_LU_decomp(V, p, &s); // faktoryzacja LU
    gsl_linalg_LU_solve(V, p, &b.vector, a); // rozwiązanie równania Ax=b

    cout << setprecision(4) << fixed; // precyzja z dokładnością do 4 cyfr po przecinku
    for (int i = 0; i < N; i++)
        cout << "a" << N - 1 - i << " = " << gsl_vector_get(a, i) << endl;

    cout << endl;
    cout << "Wielomian interpolacyjny ma postać:" << endl;
    cout << "w(x) = ";
    for (int i = 0; i < N; i++)
    {
        if (i != N - 1)
            cout << gsl_vector_get(a, i) << "x^" << N - 1 - i << "  +  ";
        else
            cout << gsl_vector_get(a, i) << "x^" << N - 1 - i << endl;
    }

    gsl_matrix_free(V);
    gsl_vector_free(a);

    return 0;
}
