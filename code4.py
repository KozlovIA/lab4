import matplotlib.pyplot as plt
import numpy as np
from math import e, pi, sqrt
from array import *


class lab04():
    def __init__(self, n = 10, t0 = 0, T = 2, y0 = 1, h=0):
        self.n = n+1
        self.h = (T-t0)/10
        self.T = T
        self.t0 = t0
        self.y0 = y0
        self.t = [0]*self.n
        for i in range(self.n):
            self.t[i] = self.t0 + self.h*i

    def output(self, y, method, error):        # Вывод данных
        print("Метод: ", method, "\tПогрешность метода: ", error, "\n№\tt\tРеальное значение\ty")
        for i in range(len(y)):
            print(i, "\t%0.1f" % self.t[i], "\t%0.16f" % self.realf(self.t[i]), "\t%0.16f" % y[i], sep="")
        print("")

    def f(self, t, y):      # Функция для расчетов
        return (t/(1+t**2))*y
    def realf(self, t):        # Функция реальных значений
        return sqrt(1+t**2)

    def explicit_Euler(self, output=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            y[i+1] = y[i] + self.h * self.f(self.t[i], y[i])
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Явный метод Эйлера", E)
        return y
        
    def Euler_Cauchy(self, output=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            y_1 = y[i] + self.h*self.f(self.t[i], y[i])
            y[i+1] = y[i] + self.h * (self.f(self.t[i], y[i]) + self.f(self.t[i+1], y_1))/2
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Метод Эйлера-Коши", E)
        return y




def main():
    tasks = lab04()
    tasks.explicit_Euler()
    tasks.Euler_Cauchy()

    #Значения явным методом Эйлера расчитаны неверно
    #A = [1]*10
    #print(A)
    #plt.plot(t)
    #plt.show()


if __name__ == "__main__":
    main()