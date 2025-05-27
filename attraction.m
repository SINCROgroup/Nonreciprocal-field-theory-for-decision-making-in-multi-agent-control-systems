function f=attraction(right_norm,x,y,dist,kh,gamma)

%Computes the attraction exerted on the herders by the targets.
%x and y are the x and y components of the distances from each herder to
%every target (possibly considering the due to delta)

% every contribution is weighted by right_norm, a matrix stating which
% herder can chase which target, which is already multiplied by a factor
% exp(-gamma |H|)

%dist is the distance from the origin of each target used to compute the
%weight of the selection rule
    weight=exp(gamma.*dist);

%the overall force is computed as a weighted average (kh regulates the strength of the interaction)

    fx=-(kh*(right_norm.*x)*weight)./((right_norm*weight));
    fy=-(kh*(right_norm.*y)*weight)./((right_norm*weight));
    f=fillmissing([fx,fy],"constant",0);

end