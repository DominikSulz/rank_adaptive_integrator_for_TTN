function [Y1] = RK_1(Y0,tau,func,t0,t1)
% This is a explicit Runge Kutta method of order 1 (i.e. the explicit Euler
% method) for matrices.

h = t1 - t0;

% eval_f = func(t0,Y0,tau);
eval_f = func(t0,Y0);
Y1 = Y0 + h*eval_f;

end