function dydt = pid_ctrl(t,y,er,er1,er2)

% T*yout*yout' + yout = Ku
% yout' = dyout/dt
% u = Kp*er1 + Ki*er + Kd* er2
% y = [yout yout]

T=184.8891; K=77.8135; tau= 3.691;   % System parametrs
%Kp = 0.5; Ki = 0.05; Kd = 3; 
% PID Controller parametrs by Cohen-Coon method
Kp = (1.35*(tau/t)^(-1) + 0.27)/K; Ki = Kp/(((2.5*tau/T+ 0.5*(tau/T)^2)/(1 + 0.6*tau/T))*T); Kd = Kp* (0.37*(tau/T)/(1 + 0.2*(tau/T)))*T; 

dydt = [y(2); ((1/T)*(-y(2)+ K*Kd*er2 + K*Kp*er1 + K*Ki*er))];