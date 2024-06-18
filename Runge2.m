function Runge2

global A b gamma m

h = 0.01;
A = 1.3;
b = 0.09441;
gamma = 0.027;
m = 0.68;

s = 0:h:100000;
d = 0:h:100000;
f = 0:h:100000;
N = length(s);

function [ OUT ] = FUN( Y )

S=Y(1);
D=Y(2);
f=Y(3);

OUT(1) = -S*[1+m*cos(f)-D];
OUT(2) = -gamma*[D-A+S*D];
OUT(3) = b;

end

Y=zeros(1,3);
Y(1,1)=0.01;
Y(1,2)=A;
Y(1,3)=0;

Y1=Y(1,1);
Y2=Y(1,2);
Y3=Y(1,3);
t=0;

for j=1:N

   k1=FUN(Y);
   k2=FUN(Y+0.5*h*k1);
   k3=FUN(Y+0.5*h*k2);
   k4=FUN(Y+h*k3);
   Y=Y+h/6*(k1+2*k2+2*k3+k4);

   Y1(j)=Y(1,1);
   Y2(j)=Y(1,2);
   Y3(j)=Y(1,3);
   t(j)=j*h;

end
figure; plot(Y2,Y1,'k'); xlabel('D'); ylabel('S')
end
