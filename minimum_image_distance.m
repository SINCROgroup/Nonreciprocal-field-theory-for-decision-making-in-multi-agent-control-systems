function d=minimum_image_distance(x,y,L,correction)

% if correction==1 (periodic domain), corrects the distances based on the
% minimum image convention

if correction==1

        d=x-y-L*round((x-y)./(L));

else
    d=x-y;
end
end