function M = wrapMat(X,s)

M = zeros(s);
c = 1;

for i = 1:s
    for j = i+1:s
        M(i,j) = X(c);
        c = c+1;
    end
end

M = M+M';
end
