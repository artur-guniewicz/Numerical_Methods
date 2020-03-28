#
#************************************************************
#*                                                          *
#*               Artur Guniewicz - Zadanie 17               *
#*                                                          *
#************************************************************
#


# minimalizacja metodą gradientow sprzezonych

import numpy                                # operacje macierzowe
from scipy.optimize import minimize         # metoda szukania minimum
import random                               # losowanie
import matplotlib.pyplot as plt             # wykres
from mpl_toolkits.mplot3d import Axes3D     # wykres 3d
from matplotlib import cm                   # kolorowanie wykresu
import sys                                  # argumenty wywolania programu


# zwraca wartosc funkcji rosenbrocka
def rosenbrock(p):
    return((1-p[0])*(1-p[0]) + 100*(p[1]-p[0]*p[0])*(p[1]-p[0]*p[0]))


# dodaje punkt do listy krokow metody CG
def add_point(x):
    global steps
    steps.append(x)


# globalna zmienna pokazujaca kolejne kroki minimalizacji
steps=[]

if __name__ == "__main__":
    if(len(sys.argv)!=4):
        print("Opcje wywolania programu")
        print("<nazwa> <liczba punktow do wylosowania> <granica> <tryb>")
        print("<granica>- w jakich granicach program ma rysowac wykres- na przyklad 20 da granice [-20,20]")
        print("<tryb> calc-tylko wylicza minima bez szkicowania wykresu")
        print("<tryb> show-wylicza minima i szkicuje droge minimalizacji(tylko do 10 punktow)")


    # 10 punktow maksymalnie dla trybu rysowania
    ile_punktow=int(sys.argv[1])
    border=int(sys.argv[2])

    # przygotowanie do wykresu
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # rysowanie funkcji rosenbrocka
    if sys.argv[3]=="show":
        X = numpy.linspace(-border,border,1000)
        Y = numpy.linspace(-border,border,1000)
        X, Y = numpy.meshgrid(X, Y)
        Z = rosenbrock([X,Y])
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

color=['black','green','gray','red','cyan','magenta','yellow','white','blue','o range']

    for i in range(ile_punktow):
        # wylosowanie punktu
        p=[random.uniform(-border,border),random.uniform(-border,border)]
        print(f"Wylosowano punkt ({p[0]},{p[1]})")
        steps.append(p)
       
        # zastosowanie minimalizacji metoda gradientow sprzezonych
        found=minimize(rosenbrock,p,method='CG',callback=add_point)
        print(f"Znalezione minimum ({found.x[0]}, {found.x[1]}) wartosc funkcji wynosi= {rosenbrock(found.x)}, liczba iteracji={len(steps)}")
        print()
        
        # dodawanie sciezki do wykresu
        if(sys.argv[3]=="show"):
            temp_x=[i[0] for i in steps]
            temp_y=[i[1] for i in steps]
            temp_z=[rosenbrock([temp_x[i],temp_y[i]]) for i in range(0,len(temp_x))]
            ax.plot3D(temp_x,temp_y,temp_z,color=color[i], linestyle='dashed', linewidth=2, markersize=12)
        
        steps.clear()

    # wyswietlenie wykresu 3D
    if sys.argv[3]=="show":
        plt.show()
