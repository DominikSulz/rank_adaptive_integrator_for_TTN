function [Y,tau] = rand_TTN_complex(n_l,r_l,r_tau)
% This function creates a random TTN with wanted sizes

%% Y
% generate the TTN
Y = cell(1,5);
tau = cell(1,3);
r_full = [r_tau,1];
C0 = tensor(rand(r_full),r_full) + 1i*tensor(rand(r_full),r_full);
Mat = tenmat(C0,4,1:3);
[Mat,~] = qr(double(Mat).',0);
C0 = tensor(mat2tens(Mat.',r_tau,4,1:3),[r_tau 1]);
% C0 = tensor(zeros(r_full),r_full);
% for i=1:min(r_tau)
%     C0(i,i,i,1) = 10^-i;
% end
Y{end-1} = 1; %identity
Y{end} = C0;

% tau1
r_tau1 = [r_l(1),r_l(3),r_l(5),r_tau(1)];
tau1 = cell(1,length(r_tau1)+1);

C0 = tensor(rand(r_tau1)) + 1i*tensor(rand(r_tau1));
Mat = tenmat(C0,4,1:3);
[Mat,~] = qr(double(Mat).',0);
C0 = mat2tens(Mat.',r_tau1,4,1:3);

tau1{end} = tensor(C0);
tau1{4} = eye(r_tau(1),r_tau(1));
tau1{1} = orth(rand(n_l(1),r_l(1)) + 1i*rand(n_l(1),r_l(1)));
tau1{2} = orth(rand(n_l(3),r_l(3)) + 1i*rand(n_l(3),r_l(3)));
tau1{3} = orth(rand(n_l(5),r_l(5)) + 1i*rand(n_l(5),r_l(5)));
Y{1} = tau1;
% tree structure
tau{1} = cell(1,3);

% tau2
r_tau2 = [r_l(4),r_l(2),r_tau(2)];
tau2 = cell(1,length(r_tau2)+1);

C0 = tensor(rand(r_tau2)) + 1i*tensor(rand(r_tau2));
Mat = tenmat(C0,3,1:2);
[Mat,~] = qr(double(Mat).',0);
C0 = mat2tens(Mat.',r_tau2,3,1:2);

tau2{end} = tensor(C0);
tau2{3} = eye(r_tau(2),r_tau(2));
tau2{1} = orth(rand(n_l(4),r_l(4)) + 1i*rand(n_l(4),r_l(4)));
tau2{2} = orth(rand(n_l(2),r_l(2)) + 1i*rand(n_l(2),r_l(2)));
Y{2} = tau2;
% tree structure
tau{2} = cell(1,2);

% tau3
Y{3} = orth(rand(n_l(6),r_l(6)) + 1i*rand(n_l(6),r_l(6)));

end

