% Nonlinear least squares optimization problem

function F = first_order_fit_function_4(x,y,t)

theta = x(1);
K = x(2);
tau = x(3);
beta = x(4);

sys_fit = tf(K*[-beta 1],[tau 1]);
sys_fit.inputdelay = theta;
[y_fit,t_fit] = step(sys_fit,t);

F = y-y_fit;  
