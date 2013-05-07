%PID Controller
clear all;
close all;
clf;

T=236.1775; K=82.9; tau= 1.0630;  % System Parameters
setpoint = 10;         % setpoint, Celsius
startpoint = -30;      % startpoint
br_MaxPidControlValue = abs(60100*10^(-4)); %absolute current to maximum for the control (60100 * 100uA = 6.01 Amps).
uplimit = br_MaxPidControlValue/16;
br_MaxErrorValue = abs(80 - 273.15); %Maximum error values will be ~800 or 80 degK


Tf = 1000;                 % Simulation time
dt = 1;

sys=tf([K],[T 1],'ioDelay',tau);
dsys=c2d(sys,dt,'z');
[num,den]=tfdata(dsys,'v');

u_1=0;u_2=0; 
y_1=0;y_2=0;

x=[0,0,0]';
error_old=0;
u_old = 0;

%Kp=100;Ki=0.1;Kd=1000;
% PID Controller parametrs by Cohen-Coon method
Kp = (1.35*(tau/T)^(-1) + 0.27)/K ; Ki = Kp/(((2.5*tau/T+ 0.5*(tau/T)^2)/(1 + 0.6*tau/T))*T); Kd = Kp* (0.37*(tau/T)/(1 + 0.2*(tau/T)))*T; 
% PID Controller parametrs by Ziegler-Nichols method 
%Kp = 1.2*(tau/T)^(-1)/K; Ki = Kp/(T*2*(tau/T)); Kd = Kp*(T*0.5*(tau/T));

Kp = Kp * 30;
Ki = Ki * 100;
Kd = Kd * 100;

for k=1:1:Tf
    time(k)=k*dt;
    
    u(k)=Kp*x(1)+Kd*x(2)+Ki*x(3);   %PID Controller

    % Restricting the output of controller
    if  u(k) < 0
       u(k)= 0;
    end
    if uplimit < u(k) - u_old
       u(k)= u(k) + uplimit;
    end
    if br_MaxPidControlValue < u(k)
       u(k)= br_MaxPidControlValue;
    end

    %Linear model
    if time(k) <= tau
        yout(k)= startpoint;
    else
        yout(k)=-den(2)*y_1+num(2)*u_1;
    end
    error(k)=setpoint-yout(k);
    
    % Restricting error: Maximum error values will be ~800 or 80 degK
    if br_MaxErrorValue <= error
        error = br_MaxErrorValue;
    end
    if error <= -br_MaxErrorValue
        error = -br_MaxErrorValue;
    end

    %Return of PID parameters
    u_2=u_1;u_1=u(k);   
    y_2=y_1;y_1=yout(k);

    x(1)=error(k);                %Calculating P
    x(2)=(error(k)-error_old)/dt;   %Calculating D
    x(3)=x(3)+error(k)*dt;        %Calculating I
    
    % Restricting I
    if br_MaxPidControlValue <= x(3)
        x(3) = br_MaxPidControlValue;
    end
    if x(3) <= -br_MaxPidControlValue
        x(3) = -br_MaxPidControlValue;
    end

    error_old=error(k);
    u_old = u(k);
end

plot(time,yout,'b');
xlabel('time(s)');ylabel('yout');
grid on;


