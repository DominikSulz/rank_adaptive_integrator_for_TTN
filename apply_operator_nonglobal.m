function [B] = apply_operator_nonglobal(X,A,d)
% X is a TTN and A a operator in TTN form. B = A(X)

B = X;
m = length(X) - 2;
s = size(X{end});
for ii=1:m
    if 1==iscell(X{ii})
        B{ii} = apply_operator_nonglobal(X{ii},A{ii});
    else
        B{ii} = apply2leaf(X{ii},A{ii});
    end
end

tmp = mult_core(double(A{end}),double(X{end}));
s2 = size(tmp);
if s(end)==1
    B{end} = tensor(tmp,[s2 1]);
else
    B{end} = tensor(tmp);
    B{end-1} = eye(s2(end),s2(end));
end

end

function [Y1] = apply2leaf(Y,A)
% This function applies the operator A leaf Y

U = A;
dum = [];
s = size(U);
for ii=1:s(2)
    M = reshape(U(:,ii),sqrt(s(1)),sqrt(s(1)));
    dum = [dum M*Y];
end
Y1 = dum;

end

function [value] = mult_core(C1,C2)

s1 = size(C1);
s2 = size(C2);
nbDim = max(length(s1), length(s2));

% We do the product
value = reshape(C1,numel(C1),1) * reshape(C2,1,numel(C2));

% And put the numbers back into their place
value = reshape(value,[s1 s2]);
value = permute(value,reshape([nbDim+1:2*nbDim;1:nbDim],1,2*nbDim));
value = reshape(value, s1.*s2);

% Idea from:
% https://de.mathworks.com/matlabcentral/fileexchange/32578-superkron


% old code
% % "Multiplies" the cores C1 and C2
% tmp1 = size(C1);
% tmp2 = size(C2);
% d = length(tmp1);
% 
% switch d
%     
%     case 2
%         C = kron(C1,C2);
%         
%     case 3
%         C = zeros(tmp1.*tmp2);
%         for l=1:tmp1(end)
%             for m=1:tmp2(end)
%                 C(:,:,m+(l-1)*tmp2(end)) = ...
%                     kron(C1(:,:,l), C2(:,:,m));
%             end
%         end
%         
%     case 4
%         C = zeros(tmp1.*tmp2);
%         for l=1:tmp1(end)
%             for m=1:tmp2(end)
%                 for p=1:tmp1(end-1)
%                     for q=1:tmp2(end-1)
%                         C(:,:,q+(p-1)*tmp2(end-1),m+(l-1)*tmp2(end)) = ...
%                             kron(C1(:,:,p,l), C2(:,:,q,m));
%                     end
%                 end
%             end
%         end
% end
% value = C;
end