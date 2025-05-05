function [pre_comp,init] = precomp_Q_mult(Y0_pre,Y0,tau)
% Y0_pre is a cell with m tensors, which correspond to the pre-contractions
% at the highes node without the Q matrices. This functions adds now all
% the Q matrices and gives out a TTN. At the end of the cell is the object,
% which correspongs to the precomputed tensor, at the first m elements
% there is a cell.
% Further, the initial data is also computed in this function.

m = length(Y0_pre);
pre_comp = cell(1,m+1);
init = cell(1,m+1);
s = size(Y0{end});
if s(end) == 1 % checks if we are at the root tensor
    top = 1;
else
    top = 0;
end


for ii=1:m
    % compute the Q_tau_i and S_tau_i matrix of the node
    v = 1:m+1;
    v = v(v~=ii);
    Mat_C = tenmat(Y0{end},ii,v);
    [Q0_i,S0_i_T] = qr(double(Mat_C).',0);
    Q0_i = tensor(Q0_i);
    if top == 1
        tensor_pre = Y0_pre{ii};
    else
        tensor_pre = Y0_pre;
    end
    % multiply the Q matrix to the Y0_pre tensor
    if iscell(Y0{ii}) == 1
        pre_save = ttt(Q0_i,tensor_pre,2,1);  
        tmp = tensor(conj(double(pre_save)));
        % FIXME: Unklar welche Dimension in tensor_pre die 0-te etc. ist.
        % -> evtl. eine Liste mit abspeichern, welche die Nummerierung und
        % Permutationen trackt!?
        pre_save = ttt(tmp,Q0_i,m+1,2);
        [pre_comp{ii},init{ii}] = precomp_Q_mult(pre_save,Y0{ii},tau{ii});
        pre_comp{ii}{end+1} = pre_save; 
      
    end
    
    % get the initial data for the i-th subtree
    tmp = Ytau_i(tau{ii},Y0{ii},S0_i_T.'); % initial value 
    if iscell(Y0{ii}) == 1
        init{ii}{end} = tmp{end};
    else
        init{ii} = tmp;
    end
    
end


end