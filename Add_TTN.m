function [C] = Add_TTN(A,B,tau)
% This function calculates an approximation of A+B, where A and B are TTN's
% of the same size. The Addition is defined recursevly.

% initialize C
C = A;
% number of branches from the root. Second case when A and B are matrices
m = length(A) - 2;
% calculates the rank of A and B 
r_A = size(A{end});
r_B = size(B{end});

%% core
r_new = r_A + r_B;
if (r_A(end) + r_B(end))==2     % when we have singleton dimension
    r_new(end) = 1;
end
% S_C = zeros(r_new);
if length(r_A) == 2 || (length(r_A) == 3 && r_A(end) == 1) 
    S_C(1:r_A(1),1:r_A(2)) = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2)) = double(B{end});
elseif length(r_A) == 3 || (length(r_A) == 4 && r_A(end) == 1)
    S_C(1:r_A(1),1:r_A(2),1:r_A(3)) = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2),...
        r_A(3)+1:r_A(3)+r_B(3)) = double(B{end});
elseif length(r_A) == 4 || (length(r_A) == 5 && r_A(end) == 1)
    S_C(1:r_A(1),1:r_A(2),1:r_A(3),1:r_A(4)) = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2),...
        r_A(3)+1:r_A(3)+r_B(3),r_A(4)+1:r_A(4)+r_B(4)) = double(B{end});
elseif length(r_A) == 5 || (length(r_A) == 6 && r_A(end) == 1)
    S_C(1:r_A(1),1:r_A(2),1:r_A(3),1:r_A(4),1:r_A(5)) = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2),...
        r_A(3)+1:r_A(3)+r_B(3),r_A(4)+1:r_A(4)+r_B(4),...
        r_A(5)+1:r_A(5)+r_B(5)) = double(B{end});
elseif length(r_A) == 6 || (length(r_A) == 7 && r_A(end) == 1)
    S_C(1:r_A(1),1:r_A(2),1:r_A(3),1:r_A(4),1:r_A(5),1:r_A(6)) ...
        = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2),...
        r_A(3)+1:r_A(3)+r_B(3),r_A(4)+1:r_A(4)+r_B(4),...
        r_A(5)+1:r_A(5)+r_B(5),r_A(6)+1:r_A(6)+r_B(6)) = double(B{end});
elseif length(r_A) == 10 ||(length(r_A) == 11 && r_A(end) == 1)
    S_C(1:r_A(1),1:r_A(2),1:r_A(3),1:r_A(4),1:r_A(5),1:r_A(6),...
        1:r_A(7),1:r_A(8),1:r_A(9),1:r_A(10)) = double(A{end});
    S_C(r_A(1)+1:r_A(1)+r_B(1),r_A(2)+1:r_A(2)+r_B(2),...
        r_A(3)+1:r_A(3)+r_B(3),r_A(4)+1:r_A(4)+r_B(4),...
        r_A(5)+1:r_A(5)+r_B(5),r_A(6)+1:r_A(6)+r_B(6),...
        r_A(7)+1:r_A(7)+r_B(7),r_A(8)+1:r_A(8)+r_B(8),...
        r_A(9)+1:r_A(9)+r_B(9),r_A(10)+1:r_A(10)+r_B(10)) = double(B{end});
end
C{end} = tensor(S_C,r_new);

%% branches
for i=1:m
    if 1 == iscell(A{i}) % A is a TTN
        dum = Add_TTN(A{i},B{i},tau{i});
        [C{i},Ri] = orth_core(dum);
        C{end} = ttm(C{end},Ri,i); % orthonormalize the core
    else % A_i and B_i are basis matrices (leafs)
        [Qi,S,V] = svd([A{i} B{i}],'econ'); % davor ohne econ
        Ri = S*V';
        %         [Qi,Ri] = qr([A{i} B{i}],0);     % when truncating its better to use
        %         svd, as we want the maximum eigenvalues
        r_tilde = rank(Ri);
        if r_tilde ~= 1
            Qi = Qi(:,1:r_tilde);
            Ri = Ri(1:r_tilde,:);
        end
        C{i} = Qi;
        C{end} = ttm(C{end},Ri,i); 
    end
end
s = size(C{end});
C{end-1} = eye(s(end),s(end));

end

function [X,R] = orth_core(T)
% This function orthonormalizes the core of T

m = length(T) - 2;
Mat = double(tenmat(T{end},m+1,1:m));
[Mat,S,V] = svd(Mat.','econ'); % davor ohne econ
R = S*V';
%[Mat,R] = qr(Mat.',0);
s = size(T{end});
rC = rank(R);
if rC ~= 1
    Mat = Mat(:,1:rC);
    R = R(1:rC,:);
    s(end) = rC;
end
if size(Mat.',1) ~= prod(s(m+1)) || ...
   size(Mat.',2) ~= prod(s(1:m))
    1
end
X = T;
X{end} = tensor(mat2tens(Mat.',s,m+1,1:m),s);
X{end-1} = eye(rC,rC);
end