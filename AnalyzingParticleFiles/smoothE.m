% SMOOTH.M
% 
% Copie de smooth.m provenant de TCV !!!

function Y=smoothE(X,ns)


% Y=SMOOTHE(X,ns)
%
% gaussian smoothing of data over ns points
%       X(n,m)          data to smooth (smoothing is done along columns)
%        ns             smoothing width (ns points on each side)
%       Y(n,m)          smoothed data
%                                                       ThDdW 3/89

[n,m]=size(X);
if n==1, X=X'; [n,m]=size(X); end
ns=min(floor((n-2)/2),ns);      % ns may not exceed n/2

if ns>0
        w=(-ns:ns)/ns;
        w=exp(-3 * w.^2);       % ponderation vector
        sw=sum(w); 
        ssw=cumsum(w);
        w=w/sw;
        X(n+1:n+2*ns+1,:)=zeros(2*ns+1,m);
        for i=1:m               % smooth
                Y(:,i)=filter(w,1,X(:,i));
                end
        Y=Y(ns+1:ns+n,:);

        for j=1:ns              % renormalise edges
                Y(j,:)=Y(j,:)*sw/ssw(j+ns);
                Y(n-j+1,:)=Y(n-j+1,:)*sw/ssw(j+ns);
                end
else
        Y=X;
end


