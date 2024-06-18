import numpy as np
import matplotlib.pyplot as plt

EPS = 1e-3
def f(x):
    return x ** 3 - 1.6 * (x ** 2) - 2.4 * x + 0.3

def df(x):
    return 3 * x **2 - 3.2 * x - 2.4

x = np.arange(-3, 4, 0.01)
plt.plot(x, f(x))
plt.xlabel('Ось х')
plt.ylabel('Ось y')
plt.title('График кубической функции')
plt.grid(True)
plt.show()

print(f'Mетод половинного деления')
print(f"Введите левую границу интервала а")
a = float(input())
print(f"Введите правую границу интервала b")
b = float(input())

def method_pd(a, b):
    n = 0
    while abs(b - a) > 2 * EPS:
        c = (a + b) / 2
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
        n += 1
    root = (a + b) / 2
    ksi = f((a + b) / 2)
    return root, n, ksi

root, n, ksi = method_pd(a, b)
print(f'Корень = {root}, Количество итераций n = {n}, Невязка ksi = {ksi}')

print(f'Метод Ньютона')
print(f'Введите значение начального приближение х0:')
x_0 = float(input())
def Newton_method(a, b, x_0):
    n = 0
    while True:
        x_k = x_0 - f(x_0)/df(x_0)
        t = abs(x_k - x_0)
        if t < EPS:
            break
        x_0 = x_k
        n += 1
    root = x_k
    ksi = f(x_k)
    return root, n, ksi

root, n, ksi = Newton_method(a, b, x_0)
print(f'Корень = {root}, Количество итераций n = {n}, Невязка ksi = {ksi}')