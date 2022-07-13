function [Y1] = RK_2(Y0,tau,func,t0,t1)
% This is a Runge Kutta method of order 2 for time discretization
% of a ODE

h = t1 - t0;

eval1 = func(t0,Y0);
eval2 = func(t0+h,Y0+h*eval1);

Y1 = Y0 + 0.5*h*(eval1 + eval2);

end