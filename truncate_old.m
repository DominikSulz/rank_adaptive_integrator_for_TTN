function[Y1] = truncate_old(Y,tol,r_max,r_min)

Y1 = Y;
m = length(Y) - 2;


%% if we are at the root tensor - we must fulfill rank condition
tmp2 = size(Y{end});

if tmp2(end) == 1 
    rank_save = zeros(m,2);
    
    for ii=1:m
        v = 1:(m+1);
        v = v(v~=ii);
        tmp = double(tenmat(Y1{end},ii,v));
        [~,S,~] = svd(tmp,"econ"); 
        [s1,s2] = size(S);
%         rk = 0;
        ss = min(s1,s2);
%         count = 1;
        tol_S = tol*norm(S,'Fro'); %Änderung: davor gar nicht 
%         dum = 1;
        
        S_diag = diag(S);
        rk = [];
        for j=1:ss
            tmp = sqrt(sum(S_diag(j:ss).^2));
            if tmp < tol_S % zuvor tol_S/m
                rk = j-1;
                break
            end
        end
        if 1==isempty(rk)
            rk = ss;
        end
        rk = max(rk,r_min);
        rk = min(rk,r_max);
        
        rank_save(ii,1) = rk;
        rank_save(ii,2) = ss;
    end
    %% check if we fulfill the rank condition
    bool = 1;
    for jj=1:m
        tmp = rank_save(:,1);
        tmp(jj) = 1;
        if prod(tmp) < rank_save(jj)
            bool = 0;
        end
    end
    
    if bool == 0
        w = max(rank_save(:,1));
        if w<rank_save(:,2)
            rank_save = w*ones(size(rank_save(:,1)));
        else
            % tmp = min(rank_save(:,1));
            count = 1;
            cc = rank_save(:,1);
            while count < max(rank_save(:,2))
                w = max(cc);
                cc = cc(cc~=w);
                if cc < rank_save(:,2)
                    tmp = max(cc);
                    count = max(rank_save(:,2)) + 1;
                end
                count = count + 1;
            end
            rank_save = tmp*ones(size(rank_save(:,1)));
        end
    else
        rank_save = rank_save(:,1);
    end
    
    %% truncate
    for ii=1:m
        % truncate ii-dimension of Y1{end}
        v = 1:(m+1);
        v = v(v~=ii);
        tmp = double(tenmat(Y1{end},ii,v));
        s = size(Y1{end});
        [P,S,Q] = svd(tmp,"econ"); % hier evtl. svd(tmp,"econ"); einsetzen
        
%         S(rank_save(ii)+1:end,rank_save(ii)+1:end)
        
        P = P(:,1:rank_save(ii));
        S = S(1:rank_save(ii),1:rank_save(ii));
        Q = Q(:,1:rank_save(ii));
        s(ii) = rank_save(ii);
        Y1{end} = tensor(mat2tens(S*Q',s,ii),s);
        
        if 1==iscell(Y1{ii})
%             Y1{ii} = truncate(Y1{ii},tol,r_max,r_min);
            m3 = length(size(Y1{ii}{end}));
            Y1{ii}{end} = ttm(Y1{ii}{end},P.',m3);
            Y1{ii}{end-1} = eye(s(ii),s(ii));
%             Y1 = rounding(Y1,tau);
            Y1{ii} = truncate_old(Y1{ii},tol,r_max,r_min);
            
%             % re-orthonormalization
%             tmp = double(tenmat(Y1{ii}{end},m3,1:m3-1)).';
%             [Q,R] = qr(tmp,0);
%             s = size(Y1{ii}{end});
%             [~,s(end)] = size(Q);
%             Y1{ii}{end} = tensor(mat2tens(Q.',s,m3),s);
%             Y1{end} = ttm(Y1{end},R,ii);
            
        else
            Y1{ii} = Y1{ii}*P;
        end
%         Y1 = rounding(Y1,tau);
    end

else 
    %% case when we are not at the root tensor    
    for ii=1:m
        v = 1:(m+1);
        v = v(v~=ii);
        tmp = double(tenmat(Y1{end},ii,v));
        s = size(Y1{end});
        [P,S,Q] = svd(tmp,"econ"); % hier evtl. svd(tmp,"econ"); einsetzen
        [s1,s2] = size(S);
%         rk = 0;
        ss = min(s1,s2);
%         count = 1;
        tol_S = tol*norm(S,'Fro'); %Änderung: davor gar nicht 
%         dum = 1;
        
        S_diag = diag(S);
        rk = [];
        for j=1:ss
            tmp = sqrt(sum(S_diag(j:ss).^2));
            if tmp < tol_S % zuvor /m
                rk = j-1;
                break
            end
        end
        if 1==isempty(rk)
            rk = ss;
        end
        rk = max(rk,r_min);
        rk = min(rk,r_max);
        
%         S(rk+1:end,rk+1:end)

        P = P(:,1:rk);
        S = S(1:rk,1:rk);
        Q = Q(:,1:rk);
        s(ii) = rk;
        Y1{end} = tensor(mat2tens(S*Q',s,ii),s);
        
        % if Y{ii} is a cell
        if 1==iscell(Y1{ii}) 
%             Y1{ii} = truncate(Y1{ii},tol,r_max,r_min);
            m3 = length(size(Y1{ii}{end}));
            Y1{ii}{end} = ttm(Y1{ii}{end},P.',m3); 
            Y1{ii}{end-1} = eye(s(ii),s(ii));
%             Y1 = rounding(Y1,tau);
            Y1{ii} = truncate_old(Y1{ii},tol,r_max,r_min);

        else
            Y1{ii} = Y1{ii}*P; 
        end
%         Y1 = rounding(Y1,tau);
    end
end

end

