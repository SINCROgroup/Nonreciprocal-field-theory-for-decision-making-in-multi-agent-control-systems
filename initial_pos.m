function pos=initial_pos(N,x0,xf,y0,yf)
% Generates the positions of N agents, randomly and uniformly distributed in a (periodic) box of boundaries [x0,xf] x [y0,yf] 

    pos=zeros(N,2);
    
    for i=1:N
        pos(i,1)=(xf-x0)*rand()+x0;
        pos(i,2)=(yf-y0)*rand()+y0;
    end


end