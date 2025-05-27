function out=periodic(X,x0,xf)
% when plotting,takes the array of positions X, defined on a periodic box,
% and creates the out vector with the equivalent positions in the "main" box of the particles
% in the X positions 

d=size(X,2);
N=size(X,1);
out=zeros(N,d);
for i=1:d
    for k=1:N

    if X(k,d)>xf
%         out(k,d)=rem(X(k,d),(xf-x0))+x0;
        out(k,d)=x0+mod((X(k,d)-xf),xf-x0);
    elseif X(k,d)<x0
        out(k,d)=xf-mod((x0-X(k,d)),xf-x0);%rem(X(k,d),(xf-x0))+(xf-x0);
    else
        out(k,d)=X(k,d);

    end
    
    end

end