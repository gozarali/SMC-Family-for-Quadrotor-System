
clc;
clear;
close all;

%% initial values
t=0:0.005:40;
n=length(t);
dt=0.005;

x1 = zeros(1,n);     x2 = zeros(1,n);       x3 = zeros(1,n);
x4 = zeros(1,n);     x5 = zeros(1,n);       x6 = zeros(1,n);
x7 = zeros(1,n);     x8 = zeros(1,n);       x9 = zeros(1,n);
x10 = zeros(1,n);    x11 = zeros(1,n);      x12 = zeros(1,n);
u1 = zeros(1,n);     u2 = zeros(1,n);       u3 = zeros(1,n);
u4 = zeros(1,n);     ux = zeros(1,n);       uy = zeros(1,n);
Sphi = zeros(1,n);   Steta = zeros(1,n);    Ssai = zeros(1,n);
Sz = zeros(1,n);     Sx = zeros(1,n);       Sy = zeros(1,n);

x1(1)= pi/4;
x2(1)= pi/4;
x3(1)= pi/4;
x4(1)= pi/4;
x5(1)= pi/4;
x6(1)= pi/4;
x7(1)= 0.6;
x8(1)= 0.6;
x9(1)= 0.6;
x10(1)=0.6;
x11(1)=0.6;
x12(1)=0.45;

u1(1)=4;
u2(1)=4;
u3(1)=4;
u4(1)=4;
ux(1)=8;
uy(1)=8;

g=9.81;
m=0.650;
l=0.23;
jr=6*(10^-5);
Ixx=7.5*(10^-3);
Iyy=Ixx;
Izz=1.3*(10^-2);
b=3.13*(10^-5);
d=7.5*(10^-7);
epsilon =1;
a1=(Iyy-Izz)/Ixx;
a2=-jr/Ixx;
a3=(Izz-Ixx)/Iyy;
a4=-jr/Iyy;
a5=(Ixx-Iyy)/Izz;
b1=l/Ixx;
b2=l/Iyy;
b3=l/Izz;

Sphi(1)=10;
Steta(1)=10;
Ssai(1)=10;
Sz(1)=1;
Sx(1)=-.5;
Sy(1)=-3;

c=[6.8688 7.0164 4.9337 17.5326 13.6277 5.4284 -19.9480 -13.6229 -7.6071 -3.6624 -12.6959 -13.6463];
landa=[c(1) c(2) c(3) c(4) c(5) c(6)];
k=[c(7) c(8) c(9) c(10) c(11) c(12)];

for i=1:n-1

Sphi(i+1)=-x2(i)+landa(1)*(0-x1(i));
Steta(i+1)=-x4(i)+landa(2)*(0-x3(i));
Ssai(i+1)=-x6(i)+landa(3)*(0-x5(i));
Sx(i+1)=-x8(i)+landa(5)*(2-x7(i));
Sy(i+1)=-x10(i)+landa(6)*(2-x9(i));
Sz(i+1)=-x12(i)+landa(4)*(2-x11(i));
   
u1(i+1)=(m/(cos(x1(i))*cos(x3(i))))*(-k(4)*tanh(Sz(i)/epsilon)+g+landa(4)*(0-x12(i)));
u2(i+1)=(1/b1)*(-k(1)*tanh(Sphi(i)/epsilon)-a1*x4(i)*x6(i)-...
    x4(i)*a2*(-sqrt(abs( u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))...
    +sqrt(abs( u1(i)/(4*b)-( u2(i)/(2*b*l))+( u4(i)/(4*d))))- sqrt(abs( u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))...
    +sqrt(abs( u1(i)/(4*b)+( u2(i)/(2*b*l))+( u4(i)/(4*d)))))+landa(1)*(0-x2(i)));
u3(i+1)=(1/b2)*(-k(2)*tanh(Steta(i)/epsilon)-a3*x2(i)*x6(i)-...
    x2(i)*a4*(-sqrt(abs( u1(1)/(4*b)-( u3(i)/(2*b*l))-(u4(i)/(4*d))))...
    +sqrt(abs( u1(i)/(4*b)-( u2(i)/(2*b*l))+( u4(i)/(4*d))))- sqrt( abs(u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))...
    +sqrt(abs( u1(i)/(4*b)+( u2(i)/(2*b*l))+( u4(i)/(4*d))))) +landa(2)*(0-x4(i)));
u4(i+1)=(1/b3)*(-k(3)*tanh(Ssai(i)/epsilon)-a5*x4(i)*x2(i) +landa(3)*(0-x6(i)));
ux(i+1)=(m/u1(i))*(-k(5)*tanh(Sx(i)/epsilon) +landa(5)*(0-x8(i)));
uy(i+1)=(m/u1(i))*(-k(6)*tanh(Sy(i)/epsilon) +landa(6)*(0-x10(i)));


x1(i+1) =x1(i)+ dt*(x2(i));
x2(i+1)=x2(i)+dt*(x4(i)*x6(i)*a1+x4(i)*a2*(-sqrt(abs( u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))+...
    sqrt( abs(u1(i)/(4*b)-( u2(i)/(2*b*l))+( u4(i)/(4*d))))- sqrt( abs(u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))+...
    sqrt(abs( u1(i)/(4*b)+( u2(i)/(2*b*l))+( u4(i)/(4*d)))))+b1*u2(i));
x3(i+1)=x3(i)+ dt*(x4(i));
x4(i+1)=x4(i)+ dt*(x2(i)*x6(i)*a3+x2(i)*a4*(-sqrt(abs( u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))...
    +sqrt( abs(u1(i)/(4*b)-( u2(i)/(2*b*l))+( u4(i)/(4*d))))- sqrt( abs(u1(i)/(4*b)-( u3(i)/(2*b*l))-( u4(i)/(4*d))))+...
    sqrt(abs( u1(i)/(4*b)+( u2(i)/(2*b*l))+( u4(i)/(4*d)))))+b2*u3(i));
x5(i+1)=x5(i)+ dt*(x6(i));
x6(i+1)=x6(i)+dt*(x4(i)*x6(i)*a5+b3*u4(i));
x7(i+1)=x7(i) +dt*(x8(i));
x8(i+1)=x8(i)+dt*(ux(i)*(1/m)*u1(i)+.1+.05*sin(4*(i)*dt));
x9(i+1)=x9(i)+dt*(x10(i));
x10(i+1)=x10(i)+dt*(uy(i)*(1/m)*u1(i)+.1+.05*sin(4*(i)*dt));
x11(i+1)=x11(i)+ dt*(x12(i));
x12(i+1)=x12(i)+dt*(-g+(((cos(x1(i))*cos(x3(i))*u1(i))/m)));

J = u1(i)*u1(i)'+u2(i)*u2(i)'+u3(i)*u3(i)'+u4(i)*u4(i)'+x1*x1'+1e-20*x2(i)*x2(i)'+x3(i)*x3(i)'+(1e-20)*x4(i)*x4(i)'...
   +x5(i)*x5(i)' +(1e-20)*x6(i)*x6(i)'+(x7(i)-2)*(x7(i)-2)'+(1e-20)*x8(i)*x8(i)'+...
   +(x9(i)-2)*(x9(i)-2)'+ (1e-20)*x10(i)*x10(i)'+(x11(i)-2)*(x11(i)-2)'+(1e-20)*x12(i)*x12(i)';

end
 
figure(1)
plot(t,x1),xlabel('Time(sec)');ylabel('phi')

 
figure(2)
plot(t,x2),xlabel('Time(sec)');ylabel('p')

 
figure(3)
plot(t,x3),xlabel('Time(sec)');ylabel('theta')

 
figure(4)
plot(t,x4),xlabel('Time(sec)');ylabel('q')

 
figure(5)
plot(t,x5),xlabel('Time(sec)');ylabel('ksai')

 
figure(6)
plot(t,x6),xlabel('Time(sec)');ylabel('r')

 
figure(7)
plot(t,x7),xlabel('Time(sec)');ylabel('X')

 
figure(8)
plot(t,x8),xlabel('Time(sec)');ylabel('u')
 
figure(9)
plot(t,x9),xlabel('Time(sec)');ylabel('Y')

 
figure(10)
plot(t,x10),xlabel('Time(sec)');ylabel('V')

 
figure(11)
plot(t,x11),xlabel('Time(sec)');ylabel('Z')

 
figure(12)
plot(t,x12),xlabel('Time(sec)');ylabel('w')
 
figure(13)
plot(t,u1),xlabel('Time(sec)');ylabel('U1(N)')

figure(14)
plot(t,u2),xlabel('Time(sec)');ylabel('U2(N.m)')

figure(15)
plot(t,u3),xlabel('Time(sec)');ylabel('U3(N.m)')

figure(16)
plot(t,u4),xlabel('Time(sec)');ylabel('U4(N.m)')


