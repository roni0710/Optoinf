%Лабораторная работа №1, вариант №3.
n = 2000;
a = -10;
b = 10;
step = (b-a)/n;
x = a:step:(b-step/2);
size(x);

xx = 2 * x;
xxx = x + 5;
f = airy(2,x/3);
ff = airy(2,xx/3);
fff = airy(2,xxx/3);
figure(1);
plot(x, f);
title("График f(x) = Bi(x/3) на интервале х");
xlabel("x");
ylabel("f(x)");
figure(2);
plot(xx, ff);
title("График f(x) = Bi(x/3) на интервале 2х");
xlabel("x");
ylabel("f(x)");
figure(3);
plot(xxx, fff);
title("График f(x) = Bi(x/3) на интервале х+5");
xlabel("x");
ylabel("f(x)");

x1 = x(1:end/4);
f1 = airy(2,x1/3);
figure(4);
plot(x1, f1);
title("График первой четверти f(x) = Bi(x/3)");
xlabel("x");
ylabel("f(x)");

x2 = x(1+end/4:3*end/4);
f2 = airy(2,x2/3);
figure(5);
plot(x2, f2);
title("Половина графика f(x) = Bi(x/3)");
xlabel("x");
ylabel("f(x)");

domain = x > 0 & x < 7;
figure(6);
plot(x(domain), f(domain));
title("График на интервале (0,7) f(x) = Bi(x/3)");
xlabel("x");
ylabel("f(x)");

g = besselj(2, x);
g2 = besselj(2, x2);
F = f + g;
G = f .* g;
figure(7);
plot(x, f,"-;f(x);", x, g, "-.;g(x);", x(1:30:end), F(1:30:end), "p;f(x)+g(x);", "markersize", 2, x(1:30:end), G(1:30:end), "x-.;f(x)g(x);", "markersize", 2);
title("Много графиков");
xlabel("x");
ylabel("f(x)");

Ref = real(f);
figure(8);
plot(x, Ref);
title("График Ref(x)");
xlabel("x");
ylabel("f(x)");
Imf = imag(f);
figure(9);
plot(x, Imf);
title("График Imf(x)");
xlabel("x");
ylabel("f(x)");
Af = abs(f);
figure(10);
plot(x, Af);
title("График амплитуды f(x)");
xlabel("x");
ylabel("|f(x)|");
Ff = arg(f);
figure(11);
plot(x, Ff);
title("График фазы f(x)");
xlabel("x");
ylabel("arg(f(x))");

h = f.* exp(i*x);
Reh = real(h);
figure(12);
plot(x, Reh);
title("График Reh(x)");
xlabel("x");
ylabel("h(x)");
Imh = imag(h);
figure(13);
plot(x, Imh);
title("График Imh(x)");
xlabel("x");
ylabel("h(x)");
Ah = abs(h);
figure(14);
plot(x, Ah);
title("График |h(x)|");
xlabel("x");
ylabel("|h(x|)");
Fh = arg(h);
figure(15);
plot(x, Fh);
title("Фаза h(x)");
xlabel("x");
ylabel("arg(h(x))");

Fm1 = f.' * g;
Fm2 = f.' * g2;
FM1 = sqrt(Fm1 .* conj(Fm1));
FM2 = sqrt(Fm2 .* conj(Fm2));
figure(16);
imagesc(FM1), colormap gray, colorbar;
title("F(x,y)");
figure(17);
imagesc(FM2), colormap gray, colorbar;
title("F(x,y) половина");

figure(18);
y = x;
[X, Y] = meshgrid(x, y);
Gm = f .* Y + g .* X;
Gmcompl = conj(Gm);
GM = sqrt(Gmcompl .* Gm);
imagesc(GM), colormap gray, colorbar;
title("G(x)");

R = sqrt(X .^ 2 + Y .^2);
Phi = atan2(Y, X);
fr = airy(2, R/3);

H = fr .* exp(1i * 5 * Phi);
figure(19);
ReH = real(H);
imagesc(ReH), colormap gray, colorbar;
title("ReH(x)");
ImH = imag(H);
figure(20);
imagesc(ImH), colormap gray, colorbar;
title("ImH(x)");
absH = abs(H);
figure(21);
imagesc(absH), colormap gray, colorbar;
title("|H(x)|");
argH = arg(H);
figure(22);
imagesc(argH), colormap gray, colorbar;
title("arg(H(x))");

HH = H(end/2, :);
figure(23);
ReHH = real(HH);
plot(HH, ReHH);
title("ReH(x)");
xlabel("x");
ylabel("H(x)");
ImHH = imag(HH);
figure(24);
plot(HH, ImHH);
title("ImH(x)");
xlabel("x");
ylabel("ImH(x)");
absHH = abs(HH);
figure(25);
plot(HH, absHH);
title("|H(x)|");
xlabel("x");
ylabel("H(x)");
argHH = arg(HH);
figure(26);
plot(HH, argHH);
title("argH(x)");
xlabel("x");
ylabel("H(x)");


















































