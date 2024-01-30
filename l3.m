%Лабораторная работа №3. Вариант №11.
n = 500;
R = 5;
hr = R/n;
r = 0:hr:(R-hr);
P = 5;
hp = P/n;
p = 0:hp:(P-hp);
r1 = 0.8;
r2 = 0.9;
r3 = 3;
r4 = 4;
h1 = 10;
h2 = 0.5;

function result = fun(rr, ra, rb)
 result = (rr >= ra & rr <= rb) * 1;
endfunction

f = h1 * fun(r, r1, r2) + h2 * fun(r, r3, r4);
figure(1);
plot(r, abs(f));
title("|f|");
figure(2);
plot(r, arg(f));
title("arg(f)");

n = length(f);
[J,K] = meshgrid(1:2*n-1, 1:2*n-1);
alpha = round(sqrt((J-n).^2+(K-n).^2))+1;
r_alpha = (alpha-1) * hr;
A = zeros(2*n-1);
index_a = alpha <= n;
A(index_a) = f(alpha(index_a));
nabeg_fhaz = 3*atan2(K-n, J-n);
f_vih = A .* exp(1i*nabeg_fhaz);

figure(3);
Abs_f_vih = abs(f_vih);
imagesc(Abs_f_vih), colormap gray, colorbar ;
title("График амплитуды f вихревое");
figure(4);
Arg_f_vih = arg(f_vih);
imagesc(Arg_f_vih), colorbar, colormap gray ;
title("График фазы f вихревое");

[X,Y] = meshgrid(r, p);
K = ((2*pi/1i^(3))) * besselj(3, 2*pi*X.*(Y)).*X;
Fp = K * p .' * hp;
figure(5);
plot(p, abs(Fp));
title("|F|");
figure(6);
plot(p, arg(Fp));
title("arg(F)")

[R,T] = meshgrid(1:2*n-1, 1:2*n-1);
alpha = round(sqrt((R-n).^2+(T-n).^2))+1;
r_alpha = (alpha-1) * hr;
A2 = zeros(2*n-1);
index_a2 = alpha <= n;
A2(index_a2) = Fp(alpha(index_a2));
nabeg_fhaz2 = (-4) * atan2(T-n, R-n);
Fp_vih = A2 .* exp(1i*nabeg_fhaz2);

figure(7)
imagesc(abs(Fp_vih)), colorbar, colormap gray ;
title("График амплитуды F(p) вихревое");
figure(8);
imagesc(arg(Fp_vih)), colorbar, colormap gray;
title("График фазы F(p) вихревое");

