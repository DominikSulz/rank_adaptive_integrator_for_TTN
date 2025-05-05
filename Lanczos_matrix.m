function [Y1] = Lanczos_matrix(Y0,func,t0,t1,A,d,l_basis)

tol_basis = 10^-8;
h = t1-t0;

%% Arnoldi
[n,m] = size(Y0);
V = Y0(:)/norm(Y0(:));
eval = 0*V;
nrm = norm(Y0(:));

for k=1:(l_basis-1)
    
    % evaluation Av_k
    tmp = func(t0,reshape(V(:,k),[n m]),A,d);
    eval(:,k) = tmp(:);
    
    % computation of v_k+1
    V(:,k+1) = eval(:,k);
    for j=1:k
        h_jk = V(:,j)'*V(:,k+1); 
        V(:,k+1) = V(:,k+1) - h_jk*V(:,j);
    end
    % normalization
    if norm(V(:,k+1)) > tol_basis
        V(:,k+1) = V(:,k+1)/norm(V(:,k+1)); 
    else
        V = V(:,1:end-1);
        break
    end
end

[~,len_V] = size(V);
if len_V == l_basis
    tmp = func(t0,reshape(V(:,end),[n m]),A,d);
    eval(:,end+1) = tmp(:);
end

% compute inner products
[~,len] = size(V);
H = zeros(len,len);
M = zeros(len,len);

for ii=1:len
    for jj=1:len
        H(ii,jj) = V(:,ii)'*eval(:,jj);
        if ii==jj
            M(ii,ii) = V(:,ii)'*V(:,ii);
        else
            M(ii,jj) = V(:,ii)'*V(:,jj);
        end
    end
end

% solving the ODE
L = chol(M,'lower');
inverse = inv(L);
e = eye(len,len);
a1 = expm(-1i * inverse*H*inverse' * h)*(L'*e(:,1));
a1 = inverse'*a1; % back-transformation

% Linear combination of basis
Y1 = a1(1)*V(:,1);
for ii=2:len
    Y1 = Y1 + a1(ii)*V(:,ii);
end

% respape
Y1 = reshape(Y1,[n m]);
Y1 = nrm*Y1; % normalize back

end