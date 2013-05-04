%PID Controller
clear all;
close all;

T=184.8891; K=77.8135; tau= 3.691;  % System Parameters
setpoint = 2430*0.1 - 273.15;         % setpoint, Celsius
Tf = 200;                 % Simulation time
dt = 1;

sys=tf([K],[T 1]);
dsys=c2d(sys,dt,'z');
[num,den]=tfdata(dsys,'v');

u_1=0.0;u_2=0.0;u_3=0.0;
y_1=0;y_2=0;y_3=0;

x=[0,0,0]';
error_1=0;

for k=1:1:1200
time(k)=k*dt;

kp=1.0;ki=2.0;kd=0.01;
% PID Controller parametrs by Cohen-Coon method
%Kp = (1.35*(tau/t)^(-1) + 0.27)/K; Ki = Kp/(((2.5*tau/T+ 0.5*(tau/T)^2)/(1 + 0.6*tau/T))*T); Kd = Kp* (0.37*(tau/T)/(1 + 0.2*(tau/T)))*T; 
% PID Controller parametrs by Ziegler-Nichols method 
%Kp = 1.2*(tau/T)^(-1)/K; Ki = Kp/(T*2*(tau/T)); Kd = Kp*(T*0.5*(tau/T));
   
u(k)=kp*x(1)+kd*x(2)+ki*x(3);   %PID Controller

% Restricting the output of controller
if u(k)>= 60100
   u(k)= 60100;
end
if u(k)<= 0
   u(k)= 0;
end

%Linear model
yout(k)=-den(1)*y_1-den(2)*y_2+num(1)*u_1+num(2)*u_2;
error(k)=setpoint-yout(k);

u_3=u_2;u_2=u_1;u_1=u(k);
y_3=y_2;y_2=y_1;y_1=yout(k);
   
x(1)=error(k);                %Calculating P
x(2)=(error(k)-error_1)/dt;   %Calculating D
x(3)=x(3)+error(k)*dt;        %Calculating I
xi(k)=x(3);

error_1=error(k);
end
plot(time,yout,'b');
xlabel('time(s)');ylabel('yout');



