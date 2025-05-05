function [Y0_pre] = Y0_precontr(Y0,A)
% This function precomputes all contractions of the TTNO A with the initial
% data Y0. If we are at the root, we obtain a cell array with all
% precontracted subtrees, without appliying the Q_tau_i. If we are at a 
% lower level in the tree, we obtain a tensor, where the tensors including
% the core of A and the cores of Y0 are contracted.
% 
% note: C=2x3x4 and D=5x6x2. Then ttt(D,C,3,1) = 5x6x3x4 -> what I want

m = length(Y0) - 2;
s = size(Y0{end});
Y0_pre = cell(1,m);
count = 0; % count the number of branches added by the contractions

if s(end) == 1 % checks if we are at the root tensor
    top = 1;
else
    top = 0;
end
core = A{end};

for ii=1:m
    if iscell(Y0{ii}) == 0 % leaf
        %% computes contractions for leaves
        [n2,r] = size(A{ii});
        n = sqrt(n2);
        re_A = mat2tens(A{ii},[n n r],[1 2],3);
        
        % contract leaves of Y0 with leaf of A
        Y0_pre{ii} = ttm(re_A,{Y0{ii}',Y0{ii}'},1:2);
        
        if top == 0 % contract the result to the corresponding core of A
            tmp = tensor(conj(double(Y0_pre{ii})));
            core = ttt(tmp,core,3,ii+count);
            count = count + 1;
        end

    else % higher level
        %% recursion of contractions
        Y0_pre{ii} = Y0_precontr(Y0{ii},A{ii});
        
        if top == 0 % contract the result to the corresponding core of A
            tmp = tensor(conj(double(Y0_pre{ii})));
            m1 = length(Y0{ii}) - 2;
            core = ttt(tmp,core,m1+1,ii+count);
            count = count + 1;
        end
    end
    
end


if top == 0 % multiply the core of Y0 to the result
    m2 = length(size(core)) - 1;
    core = ttt(Y0{end},core,1:m,1:2:(m2-1)); % complex conjugation wanted
    tmp = tensor(conj(double(Y0{end})));
    len = length(size(core));
    core = ttt(tmp,core,1:m,2:(len-1)); % no complex conjugation wanted
    Y0_pre = core;
    
else % if we are at the root tensor
    pre_save = Y0_pre;
    % multiply the pre_save to core of A
    for ii=1:m
        core = A{end};
        for jj=1:m
            if ii ~= jj
                c1 = tensor(conj(double(pre_save{jj}))); % does this make sense?
                m1 = length(size(pre_save{jj}));
                core = ttt(c1,core,m1,jj);
            end
        end
        Y0_pre{ii} = tensor(double(core));
    end
end

end


% % old
% s = size(Y0{end});
% core = A{end};
% for ii=1:m
%     m1 = length(size(Y0_pre{ii}));
%     c1 = tensor(conj(double(core))); % does this make sense?
%     core = ttt(c1,core,m1,ii);
% end
% % multiply the 2 core tensors from Y0 to the core
% m1 = lenght(size(Y0{end}));
% core = ttt(Y0{end},core,1:(m1-1),1:2:(2*m-1));
% c1 = tensor(conj(double(core))); % does this make sense?
% core = ttt(c1,Y0{end},2:2:(2*m),1:(m1-1));
% Y0_pre{end} = core;
% 
% % contract cores of Y0{end} to Y0_pre{end}
% m1 = length(size(Y0_pre{end}));
% m2 = length(size(Y0{end}));
% Y0_pre{end} = ttt(Y0_pre{end},Y0{end},1:2:(m1-2),1:(m2-1)); % FIXME Check for complex conjugation?!
% Y0_pre{end} = ttt(Y0_pre{end},Y0{end},1:(m2-1),1:(m2-1));
% 
% 
% %% contractions where all except one subtree is contracted
% for ii=1:m
%     v = 1:m+1;
%     v = v(v~=ii);
%     Mat_C = tenmat(Y0{end},ii,v);
%     [Q0_i,S0_i_T] = qr(double(Mat_C).',0);
%     ten_Q = mat2tens(Q0_i.',size(Y0{end}),ii,v);
%     
%     if iscell(Y0{ii}) == 0
%         %% initial value
%         init{ii} = double(ttm(tensor(Y0{ii}),S0_i_T,2));
%     else
%         %% initial value
%         init{ii}{end} = ttm(Y0{end},S0_i_T,m+1);
%     end
%         
% end
