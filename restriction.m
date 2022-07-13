function [Z_res] = restriction(Z_tau,Y_tau,i,Q0_tau_i)
% This function does the restriction of a given tree \tau to a subtree
% \tau_i
% suppose that all the matricizations of tensors (U_tau, W_tau etc.) are
% not real matrices but in tucker format unless they are basis matrices
%
% Y = Y_tau_i --> input of the prolongation

m = length(Z_tau) - 2;

% calculate the small R_tau_i-matrix 
dum = cell(1,m);
for j=1:i-1
    dum{j} = Mat0Mat0(Y_tau{j},Z_tau{j});    
end
for j=i+1:m
    dum{j-1} = Mat0Mat0(Y_tau{j},Z_tau{j});
end
dum{m} = Z_tau{end-1};

v = 1:m+1;
v = v(v~=i);
Ten = ttm(Z_tau{end},dum,v);
Mat = tenmat(Ten,i,v);
R_tau_i = conj(Q0_tau_i.')*double(Mat).';

% calculate the restriction 
dum = Z_tau{i};
if 0 == iscell(dum)
    dum = double(ttm(tensor(dum),R_tau_i,2));
else
    m2 = length(dum) - 2;
    dum{end} = ttm(dum{end},R_tau_i,m2+1);
end
if 1 == iscell(dum)
    [r,~] = size(R_tau_i);
    dum{end-1} = eye(r,r);
end
Z_res = dum;

end