function [r] = max_rank(X)

m = length(X) - 2;
s = size(X{end});
r = max(s);

for ii=1:m
    if 1==iscell(X{ii})
        r_i = max_rank(X{ii});
        if r_i > r
            r = r_i;
        end
    end
end

end