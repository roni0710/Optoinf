%graphics_toolkit('gnuplot')

global A b gamma m
h = 0.01;
m = 0.68;
gama = 0.027;
beta = 0.09441;
A = 1.3;

TT=1e4;
N=ceil(TT/h);
%s = 0:h:100;
%d = 0:h:100;
%f = 0:h:100;
%N = length(s);
s = zeros(1,3);
d = zeros(1,3);
f = zeros(1,3);
%t = zeros(N,1);
S = @(s, d, f) -s * (1 + m * cos(f) - d);
D = @(s, d, f) -gama * (d - A + s * d);
F = @(s, d, f) beta;
s(1) = 0.01;
d(1) = 1.3;
f(1) = 0;
t(1) = 0;



for i = 1:(N-1)
  K1 = h * S( s(i), d(i), f(i) );
  L1 = h * D( s(i), d(i), f(i) );
  M1 = h * F( s(i), d(i), f(i) );

  K2 = h * S( s(i) + 1/2*h, d(i) + 1/2*L1 , f(i) + 1/2*M1 );
  L2 = h * D( s(i) + 1/2*K1, d(i) + 1/2*h, f(i) + 1/2*M1 );
  M2 = h * F( s(i) + 1/2*K1, d(i) + 1/2*L1, f(i) + 1/2*h );

  K3 = h * S( s(i) + 1/2*h , d(i) + 1/2*L2 ,f(i) + 1/2*M2 );
  L3 = h * D( s(i) + 1/2*K2 , d(i) + 1/2*h , f(i) + 1/2*M2 );
  M3 = h * F( s(i) + 1/2*K2 , d(i) + 1/2*L2 , f(i) + 1/2*h );

  K4 = h * S( s(i) + h, d(i) + L3, f(i) + M3 );
  L4 = h * D( s(i) + K3, d(i) + h , f(i) + M3);
  M4 = h * F( s(i) + K3, d(i) + L3, f(i) + h );

  s(i+1) = s(i) + 1/6 *( K1 + 2*K2 + 2*K3 + K4);
  d(i+1) = d(i) + 1/6 * (L1 + 2*L2 + 2*L3 + L4);
  f(i+1) = f(i) + 1/6 * (M1 + 2*M2 + 2*M3 + M4);

  t(i+1) = i * h;

end
%plot(d, s,'k')%,'Linewidth',2 )
figure; plot(d,s,'k'); xlabel('D'); ylabel('S')
