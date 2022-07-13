function [M] = Mat0Mat0_prod_real(U,V)

global a b d

M = Mat0Mat0_real(U,V);
M = (b-a)^d * M;

end
