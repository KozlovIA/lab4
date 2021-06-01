import matplotlib.pyplot as plt
import numpy as np
from math import pi
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

    def output(self, y, method):
        print("Метод: ", method, "\n№\tt\t\ty")
        for i in range(len(y)):
            print(i, "\t%0.1f" % self.t[i], "\t", y[i], sep="")
    def task_1(self):
        print(self.t)
        pass
    def r(self, t):
        return t/(1+t**2)
    def explicit_Euler(self, output=True):
        y = [self.y0]*self.n
        print(self.t)
        for i in range(self.n-1):
            y[i+1] = y[i] + self.h * self.r(self.t[i])
        if(output):
            self.output(y, "Явный метод Эйлера")
        
    def Euler_Cauchy(self):
        pass




def main():
    tasks = lab04()
    tasks.explicit_Euler()


    #A = [1]*10
    #print(A)
    #plt.plot(t)
    #plt.show()


if __name__ == "__main__":
    main()