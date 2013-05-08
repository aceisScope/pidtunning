% Nonlinear least squares optimization problem

function F = first_order_fit_function(x,y,t)

theta = x(1);
K = x(2);

sys_fit = tf([K],[1 0]);
sys_fit.inputdelay = theta;
[y_fit,t_fit] = step(sys_fit,t);

F = y-y_fit;  
