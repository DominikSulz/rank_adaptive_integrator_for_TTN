function [Y1,C0_tau_hat,C1_tau_hat] = TTN_integrator_complex_rank_adapt(tau,Y0,F_tau,t0,t1)
% This function does one time-step with the unconventional integrator for 
% TTNs in a recursuive way.
%
% Input: 
%       tau = representation of the tree
%       Y0 = TTN; initial value of the integration
%       F_tau = function of the ODE 
%       t0,t1 = t1 - t0 is the timestep-size
% Output:
%       Y1 = TTN; solution of ODE at time t1
global d

Y1 = cell(size(Y0));
m = length(Y0) - 2;
M = cell(1,m);

% parfor i=1:m kann parallelisieren, wenn keine Rekursionen vorkommen
for i=1:m
    %% subflow \Phi_i
    v = 1:m+1;
    v = v(v~=i);
    
    Mat_C = tenmat(Y0{end},i,v);
    [Q0_i,S0_i_T] = qr(double(Mat_C).',0);
    
    % K-step
    Y0_i = Ytau_i(tau{i},Y0{i},S0_i_T.'); % initial value for K-step
    
    F_tau_i = @(t,Y_tau_i) restriction(...
              F_tau(t,prolongation(Y_tau_i,Y0,i,Q0_i)),Y0,i,Q0_i);
    
    if 0 == iscell(tau{i})    % if \tau_i = l, l \in L
        Y1_i = RK_4(Y0_i,tau{i},F_tau_i,t0,t1);
    else % if \tau_i \notin L
        [Y1_i,C0_taui_hat,C1_tau_hat] = TTN_integrator_complex_rank_adapt(tau{i},Y0_i,F_tau_i,t0,t1);
    end
    
    % distiguish between leaf and TTN case
    if 1 == iscell(Y1_i)
        % orthonormalization in 0-dimension of the new core tensor
%         v = 1:length(size(Y1_i{end}))-1;
%         tmp = double(tenmat(Y1_i{end},length(v)+1,v)).';
%         [Q,~] = qr(tmp,0); 
%         s = size(Y1_i{end});
%         [~,s(end)] = size(Q);
%         Y1_i{end} = tensor(mat2tens(Q.',s,length(v)+1),s);
        
        m2 = length(Y1_i) - 2;
        core1 = double(tenmat(C0_taui_hat,m2+1,1:m2)).';
        core2 = double(tenmat(Y1_i{end},m2+1,1:m2)).';
%         core2 = double(tenmat(C1_tau_hat,m2+1,1:m2)).';
        tmp = [core1 core2];
        s = size(C0_taui_hat);
        [W_taui_hat,S,~] = svd(tmp,0);
%         [W_taui_hat,S] = qr(tmp,0); 
        rr = rank(S); % [~,rr] = size(S); --> if one really wants to double the rank
        W_taui_hat = W_taui_hat(:,1:rr);
        s(end) = rr;
%         [~,s(end)] = size( W_taui_hat);
        Y1_i{end} = tensor(mat2tens(W_taui_hat.',s,m2+1),s);
%         Y1_i = truncate_1level(Y1_i,10^-8,4,4);
    else
        Ui_hat = [Y1_i Y0{i}];        
        [Y1_i,~] = qr(Ui_hat,0);   
        
        
%         Ui_hat = [Y1_i Y0{i}];
%         [Y1_i,S,~] = svd(Ui_hat,0);
%         rr = rank(S);
%         Y1_i = Y1_i(:,1:rr);
    end
    
    M{i} = Mat0Mat0(Y1_i,Y0{i});    
   
    Y1{i} = Y1_i;    
end

%% subflow \Psi

% solve the tensor ODE
C0 = ttm(Y0{end},M,1:m);
C0_tau_hat = C0;

% F_ODE = @(C0,F_tau,U1_tau,t0,tau) func_ODE(C0,F_tau,Y1(1:m),t0,tau);
F_ODE = @(C0,F_tau,U1_tau,t0) func_ODE(C0,F_tau,Y1(1:m),t0);

Y1{end-1} = eye(size(Y0{end-1}));
Y1{end} = RK_4_tensor(C0,F_ODE,Y1(1:m),F_tau,t0,t1,tau);

% neu --> macht keinen Sinn am root tensor?!
C1_tau_hat = Y1{end};

end

function [X] = func_ODE(C,F_tau,U1,t)
% function [X] = func_ODE(C,F_tau,U1,t,tau)
% This function defines the function F_tau(C(t)X U_1) X U1^T, for the
% tensor-ODE. Here C(t) is in tucker form, C0 = C X M_i, i.e. the M_i are
% matrices.
 
% argument of F_tau
m = length(U1);
s = size(C);
N = cell(1,m+2);
N{end} = C;
N{end-1} = eye(s(end),s(end));
N(1:m) = U1;

% apply F_tau
% F = F_tau(t,N,tau);
F = F_tau(t,N);

% multipl. with U1^*
dum = cell(1,m);
for i=1:m
    dum{i} = Mat0Mat0(U1{i},F{i});
end
X = ttm(F{end},dum,1:m);


end