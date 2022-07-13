function [B] = apply_operator(X,A)
% X is a TTN and A a operator in TTN form. B = A(X)
global d

B = X;
for ii=1:d
    B = apply2leaf(B,A,ii,d);
end
B = apply2core(B,A);

end

function [Y] = apply2leaf(Y,A,k,d)
% This function applies the operator A to the k-th leaf of Y

if d==2
    U = A{k};
    dum = [];
    s = size(U);
    for ii=1:s(2)
        M = reshape(U(:,ii),sqrt(s(1)),sqrt(s(1)));
        dum = [dum M*Y{k}];
    end
    Y{k} = dum;
elseif d==3 
    if k==3
        U = A{2};
        dum = [];
        s = size(U);
        for ii=1:s(2)
            M = reshape(U(:,ii),sqrt(s(1)),sqrt(s(1)));
            dum = [dum M*Y{2}];
        end
        Y{2} = dum;
    else
        Y{1} = apply2leaf(Y{1},A{1},k,2);
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
        Y{1} = apply2leaf(Y{1},A{1},k_new,d_new);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        Y{2} = apply2leaf(Y{2},A{2},k_new,d_new);
    end
    
end     
end


function [Y] = apply2core(Y,A)
% This function applies the operator A to the TTN Y recursevly from the
% root to the leafs
s = size(Y{end});
Y{end} = mult_core(double(A{end}),double(Y{end}));
s2 = size(Y{end});
if s(end) == 1
    Y{end} = tensor(Y{end},[s2 1]);
    Y{end - 1} =  1;
else
    Y{end} = tensor(Y{end});
    Y{end - 1} =  eye(s2(end),s2(end));
end

m = length(Y) - 2;
for ii=1:m
    if 1==iscell(Y{ii})
        Y{ii} = apply2core(Y{ii},A{ii});
    end
end

end


function [value] = mult_core(C1,C2)
% "Multiplies" the cores C1 and C2
tmp1 = size(C1);
tmp2 = size(C2);
d = length(tmp1);

switch d
    
    case 2
        C = kron(C1,C2);
        
    case 3
        C = zeros(tmp1.*tmp2);
        for l=1:tmp1(end)
            for m=1:tmp2(end)
                C(:,:,m+(l-1)*tmp2(end)) = ...
                    kron(C1(:,:,l), C2(:,:,m));
            end
        end
        
    case 4
        C = zeros(tmp1.*tmp2);
        for l=1:tmp1(end)
            for m=1:tmp2(end)
                for p=1:tmp1(end-1)
                    for q=1:tmp2(end-1)
                        C(:,:,q+(p-1)*tmp2(end-1),m+(l-1)*tmp2(end)) = ...
                            kron(C1(:,:,p,l), C2(:,:,q,m));
                    end
                end
            end
        end
end
value = C;
end