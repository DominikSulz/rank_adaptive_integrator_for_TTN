function[Y1] = truncate_original_rank(Y,Y0)
% This function truncates Y back to its original rank of Y0

Y1 = Y;
m = length(Y) - 2;

for ii=1:m
    v = 1:(m+1);
    v = v(v~=ii);
    tmp = double(tenmat(Y1{end},ii,v));
    s = size(Y1{end});
    [P,S,Q] = svd(tmp);
    
    s0 = size(Y0{end});
    [sP1,sP2] = size(S);
    if sP1 < s0(ii) || sP2 < s0(ii)
        rk = min(sP1,sP2);
    else
        rk = s0(ii);
    end
%     S(rk+1:end,rk+1:end)
    P = P(:,1:rk);
    S = S(1:rk,1:rk);
    Q = Q(:,1:rk);
    s(ii) = rk;
    Y1{end} = tensor(mat2tens(S*Q',s,ii),s);
    
    % if Y{ii} is a cell
    if 1==iscell(Y1{ii})
        m3 = length(size(Y1{ii}{end}));
        Y1{ii}{end} = ttm(Y1{ii}{end},P.',m3);
        Y1{ii}{end-1} = eye(s(ii),s(ii));
        Y1{ii} = truncate_original_rank(Y1{ii},Y0{ii});
    else
        Y1{ii} = Y1{ii}*P;
    end

end
end
