function [A] = make_operator_kressner(X,B,d)
% This function creates a binary tree of the same structure as X, which
% represents a linear operator of the form ....

% Structure B: B is a cell array of dimension 3xd
% First row: matrices for part without coupling
% Second row: matrices for coupling; for k-th leaf; d-th entry is 0
% Third row: matrices for coupling; for k+1-th leaf; 1-st entry is 0

A = X;

for ii=1:d    
    s = size(B{1,ii});
    I = eye(s(1),s(1));
    
    U = [I(:) B{3,ii}(:) B{2,ii}(:) B{1,ii}(:)];
%     U = [I(:) I(:) I(:) I(:)];
    A = set_operator(A,U,ii,d);
end
% C = zeros(4,4);
% C(2,1) = 1;
% C(1,3) = 1;
% C(2,3) = 1;
% C(4,1) = 1;
% C(1,4) = 1;

% C = eye(4,4);

C = zeros(4, 4);
C(4, 1) = 1; 
C(3, 2) = 1; 
C(1, 4) = 1;
C(2, 3) = 1;
A{end} = tensor(C,[4 4 1]);
A{end-1} = 1;

end

function [A] = set_operator(Y,S,k,d)
% S is the matrix that shall be on k-th leaf

% conecting tensor 
B1 = zeros(4,4,4);
B1(1, 1, 1) = 1;
B1(2, 1, 2) = 1;
B1(1, 3, 3) = 1;
B1(4, 1, 4) = 1;
B1(3, 2, 4) = 1; 
B1(1, 4, 4) = 1;

% conecting tensor if the right children is a leaf
% B2 = zeros(4, 4, 4);
% B2(1, 1, 1) = 1;
% B2(2, 1, 2) = 1;
% B2(1, 2, 3) = 1;
% B2(4, 1, 4) = 1; 
% B2(3, 2, 4) = 1;
% B2 = B1;

A = Y;
if d==2
    A{k} = S;
    A{end-1} = eye(4,4);
    A{end} = tensor(B1);
elseif d==3 
    if k==3
        A{2} = S;
    else
        A{1} = set_operator(Y{1},S,k,2);
    end
    A{end} = tensor(B1);
    A{end-1} = eye(4,4);
else
    N = 0;
    c = d - 2^N;
    while 2^N <= c
        N = N + 1;
        c = d - 2^N;
    end
    
    if c<= (2^(N-1))
        dim_left = 2^(N-1) + c;
        dim_right = 2^(N-1);
    else
        dim_left = 2^N;
        dim_right = c;
    end

    if k <= dim_left
        k_new = k;
        d_new = dim_left;
        A{1} = set_operator(Y{1},S,k_new,d_new);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        A{2} = set_operator(Y{2},S,k_new,d_new);
    end
    A{end-1} = eye(4,4);
    A{end} = tensor(B1);
end     
end