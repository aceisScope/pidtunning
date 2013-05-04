% G(s) = (K*e(-tau*s))/(T*s+1)

clear all; 
clc;

T=184.8891; K=77.8135; tau= 3.691;  % System Parameters
setpoint = 2430*0.1 - 273.15;         % setpoint, Celsius
Tf = 300;                 % Simulation time
dt = 1;

yold = 0;yold1=0;
yp = []; ys_prime = [];
er = setpoint;             % Error (Initial error = setpoint) 
sum = er;                  % Sum of error
er1 = 0;                   % First derivative of error
er2 = 0;                   % Second derivative of error
erd = 0;
eold = setpoint; eold2 = 0;

for i=1:dt:Tf
    dtspan = [i i+dt];
    %Increment PID Control
    eold2 = eold ;
    eold = er;
    er = setpoint - yold; 
    er2 = (er - eold) + (eold2 - eold);
    er1 = er - eold; 
    
    %Discrete PID Control
    sum = sum + er*dt;
    erd = (er - eold)/ dt;
    init_con = [yold ; (yold-yold1)]; % Initial conditions for the diffirential equations
    yold1 = yold;
    [t,y]  = ode45(@pid_ctrl,dtspan,init_con,[],er,er1,er2);   %increment
    %[t,y]  = ode45(@pid_ctrl,dtspan,init_con,[],er,sum,erd);
    ys = y(length(y),1);
    if i <= tau
        ys_prime = [0 ys_prime];
    else
        ys_prime = [ys ys_prime];
    end
    yold = ys_prime(1);
    yp = [yp yold];
end

plot(yp);
xlabel('Time');
ylabel('Output');
title('Output of the system with PID Controller');
grid on;