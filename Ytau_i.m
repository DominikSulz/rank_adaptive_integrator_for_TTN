function [Y] = Ytau_i(tau,Y0,S_tau_i_0T)
% This function calculates the initial data Y0_tau_i of the
% K-step (the subflow \phi_i)

m = length(Y0) - 2;
if 0 == iscell(tau) % if \tau_i = l, i.e. a leaf
    Y = double(ttm(tensor(Y0),S_tau_i_0T.',2));
else 
    Y = Y0;
    Y{end} = ttm(Y0{end},S_tau_i_0T.',m+1);
end

end

