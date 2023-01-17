function [Y] = rounding(A,tau)
Y = A;
m = length(A) - 2;
for i=1:m
    if 1 == iscell(A{i}) % A{i} is a TTN
        dum = rounding(A{i},tau{i});
        [Y{i},Ri] = orth_core(dum);        
        Y{end} = ttm(Y{end},Ri,i); % orthonormalize the core
    else % if we are at a basis matrix
        [Qi,S,V] = svd(A{i},"econ"); % davor ohne econ
        Ri = S*V';     
%         [Qi,Ri] = qr(A{i},0);     % when truncating its better to use
%         svd, as we want the maximum eigenvalues
        r_tilde = rank(Ri);
        if r_tilde ~= 1
            Qi = Qi(:,1:r_tilde);
            Ri = Ri(1:r_tilde,:);
        end
        Y{i} = Qi;
        Y{end} = ttm(Y{end},Ri,i); 
    end
end
s = size(Y{end});
Y{end-1} = eye(s(end),s(end));

end

function [X,R] = orth_core(T)
% This function orthonormalizes the core of T

m = length(T) - 2;
Mat = double(tenmat(T{end},m+1,1:m));
[Mat,S,V] = svd(Mat.',"econ"); % davor ohne econ
R = S*V';
% [Mat,R] = qr(Mat.',0); % when truncating its better to use svd
s = size(T{end});
rC = rank(R);
if rC ~= 1
    Mat = Mat(:,1:rC);
    R = R(1:rC,:);
    s(end) = rC;
end
if size(Mat.',1) ~= prod(s(m+1)) || ...
   size(Mat.',2) ~= prod(s(1:m))
%     1
end
X = T;
X{end} = tensor(mat2tens(Mat.',s,m+1,1:m),s);
X{end-1} = eye(rC,rC);
end