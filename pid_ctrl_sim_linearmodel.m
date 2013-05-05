%PID Controller
clear all;
close all;
clf;

T=184.8891; K=77.8135; tau= 3.691;  % System Parameters
setpoint = 2430*0.1 - 273.15;         % setpoint, Celsius
Tf = 200;                 % Simulation time
dt = 1;

sys=tf([K],[T 1],'ioDelay',tau);
dsys=c2d(sys,dt,'z');
[num,den]=tfdata(dsys,'v');

u_1=0;u_2=0;
y_1=0;y_2=0;

x=[0,0,0]';
error_1=0;

for k=1:1:Tf
time(k)=k*dt;

%Kp=1.0;Ki=2.0;Kd=0.01;
% PID Controller parametrs by Cohen-Coon method
%Kp = (1.35*(tau/T)^(-1) + 0.27)/K; Ki = Kp/(((2.5*tau/T+ 0.5*(tau/T)^2)/(1 + 0.6*tau/T))*T); Kd = Kp* (0.37*(tau/T)/(1 + 0.2*(tau/T)))*T; 
% PID Controller parametrs by Ziegler-Nichols method 
Kp = 1.2*(tau/T)^(-1)/K; Ki = Kp/(T*2*(tau/T)); Kd = Kp*(T*0.5*(tau/T));
   
u(k)=Kp*x(1)+Kd*x(2)+Ki*x(3);   %PID Controller

% Restricting the output of controller
% if u(k)>= 60100
%    u(k)= 60100;
% end
% if u(k)<= 0
%    u(k)= 0;
% end

%Linear model
yout(k)=-den(2)*y_1+num(2)*u_1;
error(k)=setpoint-yout(k);

%Return of PID parameters
u_2=u_1;u_1=u(k);   
y_2=y_1;y_1=yout(k);
   
x(1)=error(k);                %Calculating P
x(2)=(error(k)-error_1)/dt;   %Calculating D
x(3)=x(3)+error(k)*dt;        %Calculating I

error_1=error(k);
end

plot(time,yout,'b');
xlabel('time(s)');ylabel('yout');
grid on;


