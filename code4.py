import matplotlib.pyplot as plt
import numpy as np
from math import e, exp, inf, pi, sqrt
import sys


class lab04():
    def __init__(self, n = 10, t0 = 0, T = 2, y0 = 1, h=0, index_h=1):
        self.n = n+1
        self.h = index_h*(T-t0)/n
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

    def output_2(self, y, method, Runge_rule=0, h=0):        # Вывод данных
        if(Runge_rule != 0):
            if h != 0:
                print("Метод: ", method, "Погрешность метода: ", Runge_rule, "Шаг h ", h) 
            else:
                R = [0]*len(Runge_rule)
                for i in range(len(Runge_rule)):
                    R[i] = abs(Runge_rule[i])
                print("Метод: ", method, "Погрешность метода: ") 
                print("\n№\tt\ty\t\tПогрешность по правилу Рунге")
                for i in range(len(y)):
                    print(i, "\t%0.2f" % self.t[i], "\t%0.2f" % y[i], "\t\t", Runge_rule[i], sep="")
        else:
            print("Метод: ", method, "\n№\tt\ty")
            for i in range(len(y)):
                print(i, "\t%0.2f" % self.t[i], "\t%0.2f" % y[i], sep="")
        print("")

    def plot(self, x, y, color, method, lineStyle='', markerScale = 5):
        plt.plot(x, y, color=color, label = method, marker='h', ms = markerScale, ls=lineStyle)

    def f(self, t, y):      # Функция для расчетов
        return (t/(1+t**2))*y

    def realf(self, t, plot=False):        # Функция реальных значений
        if(plot):
            tx = np.linspace(self.t0, self.T, 100)
            y = [0]*100
            for i in range(100):
                y[i] = sqrt(1+tx[i]**2)
            self.plot(tx, y, "blue", "Реальное значение", '-', 0)
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
            self.plot(self.t, y, "violet", "Метод Эйлера-Коши", markerScale=10)
        return y

    def improved_Euler(self, output=True, plot=True, h=0):
        y = [self.y0]*self.n; E = [0]*self.n
        if h==0: h = self.h
        for i in range(self.n-1):
            y_1_2 = y[i] + self.f(self.t[i], y[i])*h/2
            y[i+1] = y[i] + h * self.f(self.t[i]+h/2, y_1_2)
            E[i+1] = abs(self.realf(self.t[i]) - y[i])
        E = max(E)
        if(output):
            self.output(y, 'Усовершенствованный метод Эйлера', E)
        if(plot):
            self.plot(self.t, y, "cyan", "Усовершенствованный метод Эйлера", markerScale=8)
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

    def f_2(self, t, y):
        return -35*y-5.5*exp(-2*t)-3*t
    
    def Runge_Kutta_p3(self, output=True, plot=True, Runge_rule=False, eps=0):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            k1 = self.f_2(self.t[i], y[i])
            k2 = self.f_2(self.t[i]+self.h/2, y[i]+self.h*k1/2)
            k3 = self.f_2(self.t[i] + self.h, y[i] - self.h*k1 + 2*self.h*k2)
            y[i+1] = y[i] + self.h/6 * (k1 + 4*k2 + k3)
        #if(output):
        if(Runge_rule):
            if eps > 0:
                while True:
                    R = self.Runge_rule(y)
                    for i in range(len(R)):
                        R[i] = abs(R[i])
                    if max(R) <= eps:
                        print("R =", max(R))
                        if(output):
                            self.output_2(y, "Метод Рунге-Кутты 3-го порядка", max(R), self.h)
                        break
                    else:
                        print(max(R), self.h)
                        self.__init__(n=self.n*10, t0=self.t0, T=self.T, y0=self.y0)
                        y = self.Runge_Kutta_p3(output=False, plot=False, Runge_rule=True, eps=0)
                plot = True
            else:
                if(output):
                    self.output_2(y, "Метод Рунге-Кутты 3-го порядка", self.Runge_rule(y, self.h))
        else:
            if(output):
                self.output_2(y, "Метод Рунге-Кутты 3-го порядка")
        if(plot):
            self.plot(self.t, y, "red", "Метод Рунге-Кутты 3-го порядка")
        return y

    def Runge_Kutta_p4_1(self, output=True, plot=True, Runge_rule=False, eps=0):
        y = [self.y0]*self.n; E = [0]*self.n
        for i in range(self.n-1):
            k1 = self.f_2(self.t[i], y[i])
            k2 = self.f_2(self.t[i] + self.h/2, y[i] + (self.h/2)*k1)
            k3 = self.f_2(self.t[i] + self.h/2, y[i] + (self.h/2)*k2)
            k4 = self.f_2(self.t[i] + self.h, y[i] + self.h*k3)
            y[i+1] = y[i] + self.h/6 * (k1 + 2*k2 + 2*k3 + k4)
        if(Runge_rule):
            if eps > 0:
                while True:
                    y = [self.y0]*self.n
                    for i in range(0, self.n-1):
                        k1 = self.f_2(self.t[i], y[i])
                        k2 = self.f_2(self.t[i] + self.h/2, y[i] + (self.h/2)*k1)
                        k3 = self.f_2(self.t[i] + self.h/2, y[i] + (self.h/2)*k2)
                        k4 = self.f_2(self.t[i] + self.h, y[i] + self.h*k3)
                        y[i+1] = y[i] + self.h/6 * (k1 + 2*k2 + 2*k3 + k4)
                    R = self.Runge_rule(y, self.h)
                    for i in range(len(R)):
                        R[i] = abs(R[i])
                    if max(R) <= eps or self.n > 10**7:
                        if(output):
                            self.output_2(y, "Метод Рунге-Кутты 4-го порядка", max(R), self.h)
                        break
                    else:
                        print(max(R), self.h)
                        self.__init__(n=(self.n-1)*10, t0=self.t0, T=self.T, y0=self.y0)
            else:
                if(output):
                    self.output_2(y, "Метод Рунге-Кутты 4-го порядка", self.Runge_rule(y, self.h))
        else:
            if(output):
                self.output_2(y, "Метод Рунге-Кутты 4-го порядка")
        if(plot):
            self.plot(self.t, y, "red", "Метод Рунге-Кутты 4-го порядка")
        return y

    def Runge_rule(self, y, h, p=3):
        self.__init__(n=self.n-1, index_h=2, t0=self.t0, T=self.T, y0=self.y0)
        if(p==3):
            y_2h = self.Runge_Kutta_p3(False, False, False)
        if(p==2):
            y_2h = self.improved_Euler(False, False, h = h)
        Runge_err = [0]*len(y_2h)
        #print(len(Runge_err), len(y_2h), len(y))
        for i in range(int(len(y)/2)):
            Runge_err[i] = (y_2h[i]-y[i*2])/(2**p-1)
        self.__init__(n=self.n-1, index_h=1, t0=self.t0, T=self.T, y0=self.y0)
        return Runge_err

    def Runge_Kutta_p3_1(self, output=True, plot=True, Runge_rule=False, eps=0):
        y = [self.y0]*self.n
        if(eps <=0):
            for i in range(self.n-1):
                k1 = self.f_2(self.t[i], y[i])
                k2 = self.f_2(self.t[i]+self.h/2, y[i]+self.h*k1/2)
                k3 = self.f_2(self.t[i] + self.h, y[i] - self.h*k1 + 2*self.h*k2)
                y[i+1] = y[i] + self.h/6 * (k1 + 4*k2 + k3)
        #if(output):
        if(Runge_rule):
            if eps > 0:
                while True:
                    y = [self.y0]*self.n
                    for i in range(0, self.n-1):
                        k1 = self.f_2(self.t[i], y[i])
                        k2 = self.f_2(self.t[i]+self.h/2, y[i]+self.h*k1/2)
                        k3 = self.f_2(self.t[i] + self.h, y[i] - self.h*k1 + 2*self.h*k2)
                        y[i+1] = y[i] + self.h/6 * (k1 + 4*k2 + k3)
                    R = self.Runge_rule(y, self.h)
                    for i in range(len(R)):
                        R[i] = abs(R[i])
                    if max(R) <= eps or self.n > 10**7:
                        if(output):
                            self.output_2(y, "Метод Рунге-Кутты 3-го порядка", max(R), self.h)
                        break
                    else:
                        print(max(R), self.h)
                        self.__init__(n=(self.n-1)*10, t0=self.t0, T=self.T, y0=self.y0)
            else:
                if(output):
                    self.output_2(y, "Метод Рунге-Кутты 3-го порядка", self.Runge_rule(y))
        else:
            if(output):
                self.output_2(y, "Метод Рунге-Кутты 3-го порядка")
        if(plot):
            self.plot(self.t, y, "red", "Метод Рунге-Кутты 3-го порядка", markerScale=3)
        return y

    def improved_Euler_eps(self, output=True, plot=True, eps=0):
        y = [self.y0]*10**8; E = [0]*self.n; t=self.t0;  i=0; fy = [0]*100; h = self.h; error = 1; first = True; k = 0
        while t <= self.T:
            error = 1
            while abs(error) > eps:
                y_1_2 = y[i] + self.f(t, y[i])*h/2
                y[i+1] = y[i] + h * self.f(t+h/2, y_1_2)
                #if first:
                h = h/2
                jt = t + h
                for j in range(2):
                    y_1_2 = y[i+j] + self.f(t, y[i+j])*h/2
                    y[i+j+1] = y[i+j] + h * self.f(t+h/2, y_1_2)
                    jt=jt+h
                err = self.Runge_rule([y[i+2]], h=h, p=2); error = err[0]
                if abs(error) > eps and first == False:
                    h = h/2
                first = False
            if 1:#k % 10000 == 0:
                print(k, "__", y[i], y[i+1], y[i+2], y[i+3], "h =", h, "t =", t)
                y = [y[i+2]]*10**8
                #print(y[0])
            #sys.exit()
            t = t + h
            i = i+2
            k = k + 1
        if(output):
            self.output(y, 'Усовершенствованный метод Эйлера', E)
        if(plot):
            self.plot(self.t, y, "cyan", "Усовершенствованный метод Эйлера", markerScale=8)
        return y

def main():
    tasks = lab04()
    #inp = input("Введите номер задания (1, 2): ")
    inp = '3'
    if(inp == '1'):
        tasks.realf(0, True)
        tasks.explicit_Euler()
        tasks.Euler_Cauchy()
        tasks.improved_Euler()
        tasks.Runge_Kutta_p4()
        plt.xlabel("X Axis"); plt.ylabel("Y, Axis")
        plt.legend()
        plt.show()
    if(inp == '2'):
        tasks.__init__(T=1.5, y0=-5)
        #tasks.Runge_Kutta_p3()
        tasks.__init__(n=100, T=1.5, y0=-5)
        tasks.Runge_Kutta_p3(plot=False, Runge_rule=True)
        tasks.__init__(n=100, T=1.5, y0=-5)
        #tasks.Runge_Kutta_p4_1(plot=True, Runge_rule=True)#, eps=10**-4)
        plt.xlabel("X Axis"); plt.ylabel("Y, Axis")
        plt.legend()
        plt.show()
    if(inp == '3'):
        tasks.__init__(n=4, T=2)
        tasks.improved_Euler_eps(False, False, eps = 10**-6)
    print("Finish")

if __name__ == "__main__":
    main()