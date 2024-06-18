import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np

def df(x, y):
    return y / x - y ** 2 * np.log(x) * (np.log(x) + 2) / x

def f_analit(x):
    return x / (x * (np.log(x)) ** 2 + 1)

def rk_step(df, x, y, h):
    k1 = h * df(x, y)
    k2 = h * df(x + h / 2, y + k1 / 2)
    k3 = h * df(x + h / 2, y + k2 / 2)
    k4 = h * df(x + h, y + k3)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

def estimate_error(df, x, y, h):
    y_h = rk_step(df, x, y, h/2)
    y_2h = rk_step(df, x, y, h)
    return abs(y_h - y_2h) / 30

# Исходные данные
x = 1
y = 1
h = 0.01
epsilon = 1e-4

# Находим оптимальный шаг интегрирования
while True:
    error = estimate_error(df, x, y, h)
    if error < epsilon:
        break
    h /= 2

print(f"Оптимальный шаг квадратурной суммы: {h}")

def rk4_solver(df, x0, y0, x_end, h):
    x_runge = [x0]
    y_runge = [y0]
    x = x0
    y = y0
    while x <= x_end:
        k1 = h * df(x, y)
        k2 = h * df(x + h / 2, y + k1 / 2)
        k3 = h * df(x + h / 2, y + k2 / 2)
        k4 = h * df(x + h, y + k3)
        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x = x + h
        x_runge.append(x)
        y_runge.append(y)
    return x_runge, y_runge

def euler_method(df, x0, y0, x_end, h):
    x_values = [x0]
    y_values = [y0]
    x = x0
    y = y0
    while x <= x_end:
        y = y + h * df(x, y)
        x = x + h
        x_values.append(x)
        y_values.append(y)
    return x_values, y_values

# Начальные условия
x0 = 1
y0 = 1
x_end = 2
h = 0.01

x_values = np.arange(x0, x_end + h, h)
x_euler, y_euler = euler_method(df, x0, y0, x_end, h)
x_runge, y_runge = rk4_solver(df, x0, y0, x_end, h)
y_exact = f_analit(x_values)
x = np.arange(1, 2, 0.1)
y = f_analit(x)

# Вывод результатов в виде таблицы
headers = ["x", "y (Аналитическое)", "y (Эйлера)", "Погрешность (Эйлера)", "y (Рунге-Кутты 4 порядка)", "Погрешность (Рунге-Кутты 4 порядка)"]
table = []
max_diff_euler = max_diff_rk4 = float('-inf')
for i in range(min(len(x_values), len(y_exact), len(y_euler), len(y_runge))):
    diff_euler = y_exact[i] - y_euler[i] if i < len(y_euler) else None
    diff_rk4 = y_exact[i] - y_runge[i] if i < len(y_runge) else None
    table.append([x_values[i], y_exact[i], y_euler[i] if i < len(y_euler) else None,
                  diff_euler, y_runge[i] if i < len(y_runge) else None, diff_rk4])

    if diff_euler is not None:
        max_diff_euler = max(max_diff_euler,abs((diff_euler)))
    if diff_rk4 is not None:
        max_diff_rk4 = max(max_diff_rk4,(diff_rk4))
    max_estimate_error = max(estimate_error(df, x, y, h))

print(tabulate(table, headers=headers))
print(f"Максимальная погрешность метода Эйлера: {max_diff_euler}")
print(f"Максимальная погрешность метода Рунге-Кутты 4 порядка: {max_diff_rk4}")
print(f"Максимальнjt отклонение: {max_estimate_error}")

# Строим график
plt.plot(x, y, label='Аналитическое решение', linestyle='-.')
plt.plot(x_euler, y_euler, label='Метод Эйлера', linestyle='-')
plt.plot(x_runge, y_runge, label='Метод Рунге-Кутты 4-го порядка', linestyle='--')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.title('График функции f(x) = x / (x * log(x)^2 + 1)')
plt.grid(True)
plt.show()