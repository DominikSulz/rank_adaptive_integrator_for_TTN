function [C1] = RK_2_tensor(C0,func,U1,F_tau,t0,t1,tau)
% This function does a RK-step for a tensor valued ODE with the classical
% Runge Kutta method.

h = t1 - t0;

eval1 = func(C0,F_tau,U1,t0);
eval2 = func(C0+h*eval1,F_tau,U1,t0+h);

C1 = C0 + 0.5*h*(eval1 + eval2);

end