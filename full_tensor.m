function [T] = full_tensor(A)
% This function calculates the full tensor of a TTN A

if 0 == iscell(A)
    T = A;
else
    
    m = length(A) - 2;
    dum = cell(1,m);
    for i=1:m
        if 1 == iscell(A{i})
            dum{i} = full_tensor(A{i});
            m2 = length(A{i}) - 2;
            dum{i} = double(tenmat(dum{i},m2+1,1:m2)).';
%             Ai = A{i};
%             m2 = length(Ai) - 2;
%             dum{i} = ttm(Ai{end},Ai(1:m2),1:m2);
%             dum{i} = double(tenmat(dum{i},m2+1,1:m2)).';
        else
            dum{i} = A{i};
        end
    end
    T = ttm(tensor(A{end}),dum,1:m);
end
end