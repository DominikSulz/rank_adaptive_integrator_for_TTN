function [A] = make_operator(X,B,tau,n)
% Input:  X TTN
%         B cell array of dimension s x d - if entry is empty we choose
%         identity
%         n vector containing size of the full tensor
% Output: A operator in TTN format

[s,d]= size(B);

A = set_cores(X,s);
z = size(A{end});
A{end} = tensor(double(A{end}),[z 1]);
A{end} = sptensor(A{end});
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
    [Q,S,P] = svd(U,0);
    rr =rank(S);
    Q = Q(:,1:rr);
    S = S(1:rr,1:rr);
    P = P(:,1:rr);
    R = S*P';
    A = set_operator(A,Q,R,ii);

end
end


function [arr] = count_leaves(X)

m = length(X) - 2;
arr = zeros(m+1,1);

for ii=1:m
    if 1==iscell(X{ii})
        tmp = count_leaves(X{ii});
        arr(ii) = tmp(end);
    else
        arr(ii) = 1;
    end
end
arr(end) = sum(arr(1:end-1));

end

function [A] = set_operator(Y,Q,R,k)
% S is the matrix that shall be on k-th leaf

A = Y;
arr = count_leaves(Y);
m = length(Y) - 2;
arr2 = zeros(m,1);
arr2(1) = arr(1);
for ii=2:m
    arr2(ii) = arr2(ii-1) + arr(ii);
end

tmp = 0;
for ii=1:m
    if (k <= arr2(ii)) && (tmp < k)
        if 1 == iscell(Y{ii})
            A{ii} = set_operator(Y{ii},Q,R,k - tmp);
        else
            A{ii} = Q;
            A{end} = ttm(A{end},R,ii);
        end
    end
    tmp = tmp + arr(ii);
end

end

function [Y] = set_cores(X,s)

Y = X; 
z = size(X{end});
if z(end) == 1
    ll = length(z) - 1;
else
    ll = length(z);
end
if ll == 2
    Y{end} = id_tensor([s s]);
elseif ll == 3
    Y{end} = id_tensor([s s s]);
elseif ll == 4
    Y{end} = id_tensor([s s s s]);
elseif ll == 5
    Y{end} = id_tensor([s s s s s]);
elseif ll == 6
    Y{end} = id_tensor([s s s s s s]);
end

Y{end-1} = eye(s,s);

if 1 == iscell(X)
    m = length(X) - 2;
    for ii=1:m
        if iscell(Y{ii}) == 1
            Y{ii} = set_cores(Y{ii},s);
        end
    end
end

end

function [C] = id_tensor(s)

C = sptensor(s);
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