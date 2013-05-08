% Nonlinear least squares optimization problem

function F = first_order_fit_function_5(x,y,t)

theta = x(1);
K = x(2);
beta = x(3);

sys_fit = tf(K*[-beta 1],[1 0]);
sys_fit.inputdelay = theta;
[y_fit,t_fit] = step(sys_fit,t);
F = y-y_fit;  
