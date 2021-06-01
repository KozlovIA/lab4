import matplotlib.pyplot as plt
import numpy as np
from math import pi
from array import *


class lab04():
    def __init__(self, n = 10, t0 = 0, T = 2, y0 = 1, h=0):
        self.n = n
        self.h = (T-t0)/10
        self.T = T
        self.t0 = t0
        self.y0 = y0
        self.t = np.linspace(self.t0, self.T, self.n)

    def output(self, y, method):
        print("Метод: ", method, "\n№\t\tt\t\t\ty")
        for i in range(len(y)):
            print(i, "\t%0.16f" % self.t[i], "\t%0.16f" % y[i], sep="")
    def task_1(self):
        print(self.t)
        pass
    def r(self, t):
        return t/(1+t**2)
    def explicit_Euler(self, output=True):
        y = [0]*self.n
        for i in range(self.n-1):
            y[i+1] = y[i] + self.h * self.r(self.t[i])
        if(output):
            self.output(y, "Явный метод Эйлера")
        
    def Euler_Cauchy(self):
        pass




def main():
    tasks = lab04()
    tasks.explicit_Euler()
    A = 1.06426; B = 0.28187484
    #print("%.1d,%.2d" % (A, B))
    #plt.plot(t)
    #plt.show()


if __name__ == "__main__":
    main()