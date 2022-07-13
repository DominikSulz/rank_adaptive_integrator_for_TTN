function [Y_pro] = prolongation(Y,Y_tau,i,Q0_tau_i)
% This function does the prolongation of a given tree \tau to a subtree
% \tau_i
%
% Y = Y_tau_i --> input of the prolongation

m = length(Y_tau) - 2;
Y_pro = cell(1,m+2);
r_tau_all = size(Y_tau{end});

for j=1:m
    if i ~= j
        Y_pro{j} = Y_tau{j};
    else
        Y_pro{j} = Y;
    end
end
Y_pro{m+1} = Y_tau{end-1};
v = 1:m+1;
v = v(v~=i);
Y_pro{end} = tensor(mat2tens(Q0_tau_i.',r_tau_all,i,v),r_tau_all);

end