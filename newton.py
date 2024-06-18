import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.optimize import fsolve


def f(x):
    return 1 - math.cos(x) / 2

def g(y):
    return math.sin(y + 1) - 1.2

def f1(x, y):
    return np.sin(y + 1) - x - 1.2

def g1(x, y):
    return 2 * y + np.cos(x) - 2

def df_dx(x, y):
    return -1

def df_dy(x, y):
    return np.cos(y + 1)

def dg_dx(x, y):
    return -np.sin(x)

def dg_dy(x, y):
    return 2

def opredelitel(x, y):
    return 1 / (-2 + (np.cos(y + 1)) * np.sin(x))

def Newton(a, b):
    xk = []
    yk = []
    xk.append(a)
    yk.append(b)
    xk.append(xk[0] - opredelitel(xk[0], yk[0]) * (2 * f1(xk[0], yk[0]) - np.cos(yk[0] + 1) * g1(xk[0], yk[0])))
    yk.append(yk[0] - opredelitel(xk[0], yk[0]) * (np.sin(xk[0]) * f1(xk[0], yk[0]) - g1(xk[0], yk[0])))
    n = 1
    while (norm((xk[n] - xk[n - 1], yk[n] - yk[n - 1]), 1) > 0.001):

        n += 1
        xk.append(xk[n - 1] - opredelitel(xk[n - 1], yk[n - 1]) * (
                    2 * f1(xk[n - 1], yk[n - 1]) - np.cos(yk[n - 1] + 1) * g1(xk[n - 1], yk[n - 1])))
        yk.append(yk[n - 1] - opredelitel(xk[n - 1], yk[n - 1]) * (
                    np.sin(xk[n - 1]) * f1(xk[n - 1], yk[n - 1]) - g1(xk[n - 1], yk[n - 1])))
        if ((norm((xk[n] - xk[n - 1], yk[n] - yk[n - 1])) > norm((xk[n - 2] - xk[n - 1], yk[n - 2] - yk[n - 1]))) and (
                norm(f1(xk[n], yk[n])) > norm(f1(xk[n - 1], yk[n - 1])))):
            Newton(xk[n], yk[n])
            return 0
    print("метод Ньютона")
    print("x=", xk[n], "y=", yk[n], "кол-во итераций", n, "невязка", norm((f1(xk[n], yk[n]), g1(xk[n], yk[n])), 1))
    return 0

def Newton_simple(a, b):
    xk = []
    yk = []
    xk.append(a)
    yk.append(b)
    xk.append(xk[0] - opredelitel(xk[0], yk[0]) * (2 * f1(xk[0], yk[0]) - np.cos(yk[0] + 1) * g1(xk[0], yk[0])))
    yk.append(yk[0] - opredelitel(xk[0], yk[0]) * (np.sin(xk[0]) * f1(xk[0], yk[0]) - g1(xk[0], yk[0])))
    n = 1
    while (norm((xk[n] - xk[n - 1], yk[n] - yk[n - 1]), 1) > 0.001):
        n += 1
        xk.append(xk[n - 1] - opredelitel(xk[0], yk[0]) * (
                    2 * f1(xk[n - 1], yk[n - 1]) - np.cos(yk[0] + 1) * g1(xk[n - 1], yk[n - 1])))
        yk.append(
            yk[n - 1] - opredelitel(xk[0], yk[0]) * (np.sin(xk[0]) * f1(xk[n - 1], yk[n - 1]) - g1(xk[n - 1], yk[n - 1])))
        if ((norm((xk[n] - xk[n - 1], yk[n] - yk[n - 1])) > norm((xk[n - 2] - xk[n - 1], yk[n - 2] - yk[n - 1]))) and (
                norm(f1(xk[n], yk[n])) > norm(f1(xk[n - 1], yk[n - 1])))):
            Newton_simple(xk[n], yk[n])
            return 0
    print("модифицированный метод Ньютона")
    print(f"x = {xk[n]} \n y = {yk[n]} \n кол-во итераций = {n} \n невязка =  {norm((f1(xk[n], yk[n]), g1(xk[n], yk[n])))} ")
    return 0


x = np.linspace(-6, 6, 500)
y = np.linspace(-6, 6, 500)
y1 = []
x1 = []
for i in range(0, 500):
    y1.append(f(x[i]))
    x1.append(g(y[i]))
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.plot(x, y1)
plt.plot(x1, y)
plt.legend(['sin(y+1)-x-1.2', '2y+cos(x)-2'], loc=2)
plt.title("Графическое решение системы нелинейных уравнений")
plt.show()
print("Метод Ньютона и его модификация")
print("введите начальное приближение  x")
a_str = input()
a = float(a_str)
print("введите начальное приближение  y")
b_str = input()
b = float(b_str)

Newton(a, b)
Newton_simple(a, b)
