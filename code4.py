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


    def task_1(self):
        print(self.t)
        pass
    def r(self, t):
        return t/(1+t**2)
    def explicit_Euler(self):
        
        y = array('d', [0, 10])
        for i in range(self.n-1):
            self.t 
            y[i+1] = y[i] + self.h * self.r(self.t[i])
            print("1")
        return y
        
    def Euler_Cauchy(self):
        pass



def main():
    tasks = lab04()
    #tasks.task_1()

    #plt.plot(t)
    #plt.show()

if __name__ == "__main__":
    main()