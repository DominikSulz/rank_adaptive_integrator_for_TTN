function [M] = Mat0Mat0(U,V)
% Input are 2 TTNs and this function calculates the product
% conj(Mat_0(U))*Mat_0(V)' in a recursive way.

if 0 == iscell(U)   % if U and V are matrices
    M = conj(U.')*V;
else % if U and V are TTNs
    m = length(U) - 2;
    dum = cell(1,m);
    for j=1:m
        dum{j} = Mat0Mat0(V{j},U{j});
    end
    
    Ten = ttm(U{end},dum,1:m);
    Mat = double(tenmat(Ten,m+1,1:m));
    core2 = double(tenmat(V{end},m+1,1:m));
    M = conj(Mat)*core2.';
end

end