%Лабораторная работа №2 Вариант №2.
n = 2000;
m = 1000;
a = -5;
b = 5;
p = -5;
q = 5;
stepx = (b-a)/n;
stepe = (q-p)/m;
x = a:stepx:(b-stepx/2);
e = p:stepe:(q-stepe/2);
f = 1 ./ (1 + x.^2);
absf = abs(f);
argf = arg(f);
figure(1);
plot(x, absf);
title("|f|");
figure(2);
plot(x, argf);
title("arg(f)");

S = stepx * sum(f);
display(S)

A = exp(-2*pi*i.*e.'*x);
figure(3);
imagesc(abs(A)), colormap gray, colorbar;
title("|A|");
figure(4);
imagesc(arg(A)), colormap gray, colorbar;
title("arg(A)");

F = A * f .' * stepx;
figure(5);
plot(e, abs(F));
title("|F|");
figure(6);
plot(e, arg(F));
title("arg(F)");

B = A';
ff = B * F * stepe;
figure(7);
plot(x, abs(ff), e, abs(pi * exp(-2*pi*abs(e))));
title("|f|");
figure(8);
plot(x, arg(ff), e, arg(pi * exp(-2*pi*abs(e))));
title("arg(f)");

figure(9);
plot(x, absf, x, abs(ff));
title("Сравнение графиков амплитуд исходной и полученной функций");
figure(10);
plot(x, argf, x, arg(ff));
title("Сравнение графиков фаз исходной и полученной функций");







