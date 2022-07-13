function [Y,tau] = rand_Tucker_complex()
% This function creates a random TTN with wanted sizes

n_l = [11 11 11];
r_l = [5 5 5];

%% Y
% generate the TTN
Y = cell(1,5);
tau = cell(1,5);
r_full = [r_l 1];
C0 = tensor(rand(r_full),r_full) + 1i*tensor(rand(r_full),r_full);
Mat = tenmat(C0,4,1:3);
[Mat,~] = qr(double(Mat).',0);
C0 = tensor(mat2tens(Mat.',r_full,4,1:3),r_full);
Y{end-1} = 1; %identity
Y{end} = C0;

Y{1} = orth(rand(n_l(1),r_l(1)) + 1i*rand(n_l(1),r_l(1)));
Y{2} = orth(rand(n_l(2),r_l(2)) + 1i*rand(n_l(2),r_l(2)));
Y{3} = orth(rand(n_l(3),r_l(3)) + 1i*rand(n_l(3),r_l(3)));

end

