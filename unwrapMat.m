function X = unwrapMat(Mat)

c = 1;
for i = 1:size(Mat,1)-1
    
    for j = i+1:size(Mat,1) 
        
        X(c) = Mat(i,j); 
        c = c+1;
    end
end
