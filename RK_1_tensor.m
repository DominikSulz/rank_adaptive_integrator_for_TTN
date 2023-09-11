function [C1] = RK_1_tensor(C0,func,U1,F_tau,t0,t1,tau,A,d)
% This function does a explicit Euler step for a tensor valued ODE. 

h = t1 - t0;

% func_eval = func(C0,F_tau,U1,t0,tau);
func_eval = func(C0,F_tau,U1,t0,A,d);
C1 = C0 + h*func_eval;

end