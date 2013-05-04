% G(s) = (K*e(-tau*s))/(T*s+1)

clear all; 
clc;

T=184.8891; K=77.8135; tau= 3.691;  % System Parameters
setpoint = 2430*0.1 - 273.15;         % setpoint, Celsius
Tf = 1200;                 % Simulation time
yold = 0;yold1=0;
yp = []; ys_prime = [];
er = setpoint;                 % Error (Initial error = setpoint)       
er1 = 0;                   % First derivative of error
er2 = 0;                   % Second derivative of error
eold = setpoint; eold2 = 0;
dt = 1;
for i=1:dt:Tf
    dtspan = [i i+dt];
    eold2 = eold ;
    eold = er;
    er = setpoint - yold; %增量式PID控制
    er2 = er + eold2 - 2*eold;
    er1 = er - eold;
    init_con = [yold ; (yold-yold1)]; % Initial conditions for the diffirential equations
    yold1 = yold;
    [t,y]  = ode45(@pid_ctrl,dtspan,init_con,[],er,er1,er2); 
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