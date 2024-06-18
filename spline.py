import matplotlib.pyplot as plt
from sympy import symbols, diff
import numpy as np

def f(x):
    return 8 * np.pi / ((x * np.pi + 12) ** (3))

def Mn_1(a, b):
    x = symbols("x")
    f = 8 * np.pi / ((x * np.pi + 12) ** (3))
    maximum = 0
    f4 = diff(f, x, 4)
    c = np.linspace(2, 3.5, 100)
    for i in range(1000):
        maximum = max(maximum, abs(f4.subs(x, c[i])))
    return maximum

def define_hcorrect(a, b, eps):
    return float((24 * eps * f(b) / Mn_1(a, b)) ** (1 / 4))

def tabulation(a, b, h):
    n = int(round((b - a) / h)) + 1
    return np.array([a + i * h for i in range(n)]), np.array(
        [f(a + i * h) for i in range(n)])

def lagrange_polynomial(x, interval):
    i = 0
    interpolated_function = 0
    while i < 4:
        j = 0
        multiplier = 1
        while j < 4:
            if i != j:
                multiplier *= (x - interval[j]) / (interval[i] - interval[j])
            j += 1
        interpolated_function += f(interval[i]) * multiplier
        i += 1
    return interpolated_function

def Newton(x, x_int, y_int_l, h):
    t = (x - x_int[0]) / h
    dy0 = y_int_l[1] - y_int_l[0]
    dy1 = y_int_l[2] - y_int_l[1]
    dy2 = y_int_l[3] - y_int_l[2]
    dy0_2 = dy1 - dy0
    dy1_2 = dy2 - dy1
    dy0_3 = dy1_2 - dy0_2
    return (
        y_int_l[0]
        + dy0 * t
        + dy0_2 * t * (t - 1) / 2
        + dy0_3 * t * (t - 1) * (t - 2) / 6)

def parabolic_spline(x, x_int_l, y_int_n, h):
    matrix = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, h, h**2, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, h, h**2, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, h, h**2],
            [1, h, h**2, -1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, h, h**2, -1, 0, 0],
            [0, 1, 2 * h, 0, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 2 * h, 0, -1, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
        ]
    )
    vector = np.array(
        [
            [y_int_n[0]],
            [y_int_n[1]],
            [y_int_n[2]],
            [y_int_n[3]],
            [0],
            [0],
            [0],
            [0],
            [0],
        ]
    )
    solved = np.linalg.solve(matrix, vector)
    return evaluate_spline(x, x_int_l, solved.ravel(), 2)


def cubic_spline(x, x_int_l, y_int_n, h):
    matrix = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, h, h**2, h**3, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, h, h**2, h**3, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, h, h**2, h**3],
            [0, 1, 2 * h, 3 * h**2, 0, -1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 2 * h, 3 * h**2, 0, -1, 0, 0],
            [0, 0, 2, 6 * h, 0, 0, -2, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 2, 6 * h, 0, 0, -2, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 6 * h],
            [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    )
    vector = np.array(
        [
            [y_int_n[0]],
            [y_int_n[1]],
            [y_int_n[1]],
            [y_int_n[2]],
            [y_int_n[2]],
            [y_int_n[3]],
            [0],
            [0],
            [0],
            [0],
            [0],
            [0],
        ]
    )
    solved = np.linalg.solve(matrix, vector)
    return evaluate_spline(x, x_int_l, solved.ravel(), 3)

def evaluate_spline(x, x_int_l, coefficients, order):
    value = 0
    if x_int_l[0] <= x <= x_int_l[3]:
        i = 0
        while i <= order:
            n = (
                0
                if x < x_int_l[1]
                else order + 1 if x_int_l[1] <= x < x_int_l[2] else 2 * (order + 1)
            )
            n += i
            value += coefficients[n] * np.power(
                x - x_int_l[int(n / (order + 1))], n % (order + 1)
            )
            i += 1
    return value


a = 2
b = 3.5
eps = 10 ** (-4)

h = define_hcorrect(a, b, eps)
x, y = tabulation(a, b, h)

number_x = 1
x_Lagrange = np.array([x[(number_x + i)] for i in range(4)])
y_Lagrange = np.array(
    [lagrange_polynomial(x_Lagrange[i], x_Lagrange) for i in range(len(x_Lagrange))]
)
x, y = tabulation(x_Lagrange[0], x_Lagrange[3], (x_Lagrange[3] - x_Lagrange[0]) / 100)
y_int_Lagrange = np.array(
    [lagrange_polynomial(x[i], x_Lagrange) for i in range(len(x))]
)

errorL = []
for i in range(len(x)):
    errorL.append(abs(y_int_Lagrange[i] - y[i]))

y_Newton = np.array(
    [Newton(x_Lagrange[i], x_Lagrange, y_Lagrange, h) for i in range(len(x_Lagrange))]
)
y_int_Newton = np.array(
    [Newton(x[i], x_Lagrange, y_Lagrange, h) for i in range(len(x))]
)

errorN = []
for i in range(len(x)):
    errorN.append(abs(y_int_Newton[i] - y[i]))

yparabolic_spline = np.array(
    [parabolic_spline(x[i], x_Lagrange, y_Newton, h) for i in range(len(x))]
)
ycubic_spline = np.array(
    [cubic_spline(x[i], x_Lagrange, y_Newton, h) for i in range(len(x))]
)

error_parabolic = []
for i in range(len(x)):
    error_parabolic.append(abs(yparabolic_spline[i] - y[i]))
error_cubic = []
for i in range(len(x)):
    error_cubic.append(abs(ycubic_spline[i] - y[i]))

plt.grid()
plt.plot(x, y, label="Исходная функция")
plt.plot(x, y_int_Lagrange, label="Полином Лагранжа")
plt.plot(x, y_int_Newton, label="Полином Ньютона")
plt.plot(x, yparabolic_spline, label="Параболический сплайн")
plt.plot(x, ycubic_spline, label="Кубический сплайн")
plt.legend()
plt.show()
plt.grid()
plt.plot(x, error_parabolic, label="Параболический сплайн")
plt.plot(x, error_cubic, label="Кубический сплайн")
plt.plot(x, errorL, label="Полином Лагранжа")
plt.plot(x, errorN, label="Полином Ньютона")
plt.legend()
plt.show()

print(define_hcorrect(a, b, eps))
