function [Y1] = Lanczos_tensor(Y0,func,U1,F_tau,t0,t1,A,d,l_basis)

tol_basis = 10^-8;
h = t1-t0;

%% Arnoldi
s = size(Y0);
m = length(s);
tmp = double(tenmat(Y0,m,1:(m-1)));
sz = size(tmp);
tmp = tmp(:);

V = tmp/norm(tmp);
eval = 0*V;
nrm = norm(tmp);

for k=1:(l_basis-1)
    
    % evaluation Av_k
    C = mat2tens(reshape(V(:,k),sz),s,m);
    tmp = func(C,F_tau,U1,t0,A,d);
    tmp = tensor(double(tmp),s);
    tmp = double(tenmat(tmp,m,1:(m-1)));
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
    C = mat2tens(reshape(V(:,end),sz),s,m);
    tmp = func(C,F_tau,U1,t0,A,d); 
    tmp = tensor(double(tmp),s);
    tmp = double(tenmat(tmp,m,1:(m-1)));
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
Y1 = mat2tens(reshape(Y1,sz),s,m);
Y1 = tensor(double(Y1),s);
Y1 = nrm*Y1; % normalize back

end