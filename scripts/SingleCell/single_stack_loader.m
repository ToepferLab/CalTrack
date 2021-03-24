function [ S1] = single_stack_loader( filename )

    im = bfopen(filename);

    j=1;
    for i=1:(size(im{1},1))
        S1(:,:,j) = im{1}{i};
        j= j+1;
    end
    
end

