function [A] = make_operator_old(X,B,tau,n)
% Input:  X TTN
%         B cell array of dimension s x d - if entry is empty we choose
%         identity
%         n vector containing size of the full tensor
% Output: A operator in TTN format

[s,d]= size(B);

A = set_cores(X,s);
A{end} = eye(s,s);
A{end} = tensor(A{end},[s s 1]);
A{end-1} = 1;

for ii=1:d
    U = [];
    for jj=1:s
        if isempty(B{jj,ii}) == 1
            tmp = eye(n(ii),n(ii));
        else
            tmp = B{jj,ii}(:);
        end
        U = [U tmp(:)];
    end
    [Q,S,P] = svd(U);
    rr =rank(S);
    Q = Q(:,1:rr);
    S = S(1:rr,1:rr);
    P = P(:,1:rr);
    R = S*P';
    A = set_operator(A,Q,R,ii,d);
end
% A = rounding(A,tau);

end

function [A] = set_operator(Y,Q,R,k,d)
% S is the matrix that shall be on k-th leaf
A = Y;
if d==2
    A{k} = Q;
    A{end} = ttm(A{end},R,k); % hier unterschied sptensor zu tensor
elseif d==3 
    if k==3
        A{2} = Q;
        A{end} = ttm(A{end},R,2);
    else
        A{1} = set_operator(Y{1},Q,R,k,2);
    end
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
        A{1} = set_operator(Y{1},Q,R,k_new,d_new);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        A{2} = set_operator(Y{2},Q,R,k_new,d_new);
    end
end     
end

function [Y] = set_cores(X,s)

Y = X; 
Y{end} = id_tensor([s s s]);
Y{end-1} = eye(s,s);
if 1 == iscell(X)
    m = length(Y) - 2;
    for ii=1:m
        if iscell(Y{ii}) == 1
            Y{ii} = set_cores(Y{ii},s);
        else
            
        end
    end
end

end

function [C] = id_tensor(s)

C = zeros(s);
C = tensor(C);

l = length(s);
for ii=1:s(1)
    if l == 2
        C(ii,ii) = 1;
    elseif l == 3
        C(ii,ii,ii) = 1;
    elseif l == 4
        C(ii,ii,ii,ii) = 1;
    elseif l == 5
        C(ii,ii,ii,ii,ii) = 1;
    elseif l == 6
        C(ii,ii,ii,ii,ii,ii) = 1;    
    end
end

end