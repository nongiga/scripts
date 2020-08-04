function h = my_imagesc(img) 
% function h = my_imagesc(img) 
% same as imagesc but returns individual handles to each pixel
%
h = [] ;
for i = 1:size(img,1) 
    for j = 1:size(img,2) 
        h(i,j) = patch(j-0.5+[0 1 1 0],i-0.5+[0 0 1 1],img(i,j)+[0 0 0 0]) ;
    end
end

end