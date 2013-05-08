% Nonlinear least squares optimization problem

function F = second_order_fit_function(x,y,t)

theta = x(1);
K = x(2);
tau_1 = x(3);
tau_2 = x(4);

sys_fit = tf([K],conv([tau_1 1],[tau_2 1]));
sys_fit.inputdelay = theta;
[y_fit,t_fit] = step(sys_fit,t);

F = y-y_fit;  