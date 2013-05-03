function yprime = pid_ctrl(t,y,er,er1,er2)

T=184.8891; K=77.8135;   % System parametrs
Kp = 0.5; Ki = 0.05; Kd = 3; % PID Controller parametrs

yprime = [y(2); ((1/T)*(-y(2)+ K*Kd*er2 + K*Kp*er1 + K*Ki*er))];