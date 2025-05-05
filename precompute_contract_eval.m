function [pre_contr,initial_values] = precompute_contract_eval(Y0,Y0_init,tau,A,core_top,Y0_pre)

% contr_top is the precomputed contraction from one level above. If it is
% empty then we are at the root tensor

m = length(Y0) - 2;
pre_contr = cell(1,m+2);
contr_Y0 = cell(1,m+2);
mat_save = cell(1,m);
initial_values = cell(1,m+2);

if isempty(core_top) == 0
    pre_contr{end} = core_top;
end

%% computes the precontraction (top bottom)
for ii=1:m
    v = 1:m+1;
    v = v(v~=ii);
    Mat_C = tenmat(Y0{end},ii,v);
    [Q0_i,S0_i_T] = qr(double(Mat_C).',0);
    
    %% precomputes all initial values for all subtrees \tau for the
    % recursion. (i.e. not for the core updates)
    if iscell(Y0) == 0
        initial_values{ii} = Ytau_i(tau{ii},Y0_init{ii},S0_i_T.');
        Y0_init{ii} = initial_values{ii};
    else
        initial_values{ii}{end} = Ytau_i(tau{ii},Y0_init{ii},S0_i_T.');
        Y0_init{ii}{end} = initial_values{ii}{end};
    end
    
    %% computes the pre_contr from btop to bottom
    if iscell(Y0{ii}) == 0  % at a intermediate tensor
        core = mat2tens(Q0_i,size(Y0{end}),ii,v);
        tmp = cell(1,m-1);
        count = 1;
        for jj=1:m
            if jj ~= ii
                if iscell(contr_Y0{jj}) == 1
                    tmp{count} = contr_Y0{jj}{end};
                else
                    tmp{count} = contr_Y0{jj};
                end
                count = count + 1;
            end
        end
        tmp{end} = double(tenmat(core_top,ii,v)).';
        pre_contr{ii}{end} = ttm(core,tmp,v);
        
        % recursion
        core_top = pre_contr{ii}{end};
        [pre_contr,initial_values] = precompute_contract_eval(Y0{ii},Y0_init{ii},tau{ii},A{ii},core_top,Y0_pre{ii});
    end
    
end

end