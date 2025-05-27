function [rho,theta]=cartesian_to_polar(x,y)

% converts bidimensional (arrays) of Cartesian coordinates [x,y], into the
% corresponding array of polar coordinates

% X must be the coordinates intended as X-O under the minimum image
% convention, where O is the origin

    rho=zeros(size(x,1),1);
    theta=zeros(size(x,1),1);
    for k=1:size(x,1)
        
        rho(k)=((x(k)^2+ y(k)^2 )^.5);

        if y(k)<0
            theta(k)=2*pi-acos(x(k)/rho(k));
        else 
           theta(k)=acos(x(k)/rho(k));
        end

    end

end

    
