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

    def plot(self, x, y, color, method, lineStyle='', markerScale = 5):
        plt.plot(self.t, y, color=color, label = method, marker='h', ms = markerScale, ls=lineStyle)
        plt.xlabel="X Axis"; plt.ylabel="Y Axis"
        

    def f(self, t, y):      # Функция для расчетов
        return (t/(1+t**2))*y
    def realf(self, t, plot=False):        # Функция реальных значений
        if(plot):
            tx = np.linspace(self.t0, self.T, 100)
            y = [0]*100
            for i in range(100):
                y[i] = sqrt(1+tx[i]**2)
            #self.plot()
        return sqrt(1+t**2)

    def explicit_Euler(self, output=True, plot=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            y[i+1] = y[i] + self.h * self.f(self.t[i], y[i])
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Явный метод Эйлера", E)
        if(plot):
            self.plot(self.t, y, "green", "Явный метод Эйлера")
        return y

    def Euler_Cauchy(self, output=True, plot=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            y_1 = y[i] + self.h*self.f(self.t[i], y[i])
            y[i+1] = y[i] + self.h * (self.f(self.t[i], y[i]) + self.f(self.t[i+1], y_1))/2
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Метод Эйлера-Коши", E)
        if(plot):
            self.plot(self.t, y, "violet", "Метод Эйлера-Коши")
        return y

    def improved_Euler(self, output=True, plot=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            y_1_2 = y[i] + self.f(self.t[i], y[i])*self.h/2
            y[i+1] = y[i] + self.h * self.f(self.t[i]+self.h/2, y_1_2)
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Усовершенствованный метод Эйлера", E)
        if(plot):
            self.plot(self.t, y, "cyan", "Усовершенствованный метод Эйлера")
        return y

    def Runge_Kutta_p4(self, output=True, plot=True):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            k1 = self.f(self.t[i], y[i])
            k2 = self.f(self.t[i] + self.h/2, y[i] + (self.h/2)*k1)
            k3 = self.f(self.t[i] + self.h/2, y[i] + (self.h/2)*k2)
            k4 = self.f(self.t[i] + self.h, y[i] + self.h*k3)
            y[i+1] = y[i] + self.h/6 * (k1 + 2*k2 + 2*k3 + k4)
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, "Метод Рунге-Кутты 4-го порядка", E)
        if(plot):
            self.plot(self.t, y, "red", "Метод Рунге-Кутты 4-го порядка")
        return y


def main():
    tasks = lab04()
    tasks.realf(0, True)
    tasks.explicit_Euler()
    tasks.Euler_Cauchy()
    tasks.improved_Euler()
    tasks.Runge_Kutta_p4()
    plt.show()
if __name__ == "__main__":
    main()