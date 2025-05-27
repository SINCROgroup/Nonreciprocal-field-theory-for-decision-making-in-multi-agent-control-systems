function f=repulsion(d,dx,dy,k,r)
% computes the interaction vector f given by a repulsion whose intensity
% decays linearly with distance d, interaction range r, and coupling
% strength k.
%dx and dy are the arrays of x and y cartesian components of the relative
%displacement between any pair of agents  (dx^2+dy^2=d^2)

    fx=k*(( (r-d)./d ).*dx)  .*(d>0 & d<=r);
    fy=k*(( (r-d)./d ).*dy)  .*(d>0 & d<=r);
    f=fillmissing(cat(3,fx,fy),"constant",0);

end