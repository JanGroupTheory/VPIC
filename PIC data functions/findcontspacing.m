function [contvals] = findcontspacing(fluxvals)

    nstep = 20;
    Nconts = 100;
    eps = 0.08;

    %test = sortrows([reshape(fluxvals(2:end-1,2:end-1),1,prod(size(fluxvals)-2));1:prod(size(fluxvals)-2)]');
    
    %sqrt((fluxvals(2:end-1,3:end)-fluxvals(2:end-1,1:end-2)).^2 + (fluxvals(3:end,2:end-1)-fluxvals(1:end-2,2:end-1)).^2);
    
    dfluxvals = sqrt( (fluxvals(2:end-1,3:end)-fluxvals(2:end-1,1:end-2)).^2 + (fluxvals(3:end,2:end-1)-fluxvals(1:end-2,2:end-1)).^2);
    
    fluxrange = [minall(fluxvals) maxall(fluxvals)];
    dflux = fluxrange(2)-fluxrange(1);
   
    fluxsteps = fluxrange(1) + (0:(1/nstep):1)*dflux;
    
    gradmin = zeros(1,nstep);
    
    for i = 1:nstep
        gradmin(i) = minall(dfluxvals((fluxvals(2:end-1,2:end-1) >= fluxsteps(i))&(fluxvals(2:end-1,2:end-1) <= fluxsteps(i+1))));
    end
    
    a = (nstep)/(2*dflux*eps);
    %a = Nconts/dflux;
    q = 1;
    
    nconts = max(round(1./(gradmin*a+eps)),1);
    
    %[q sum(nconts) a]
    
    while((abs(sum(nconts)- Nconts) > 2) && (q < 10))
        
        a = a*(sum(nconts)/Nconts)^2;
        nconts = max(round(1./(gradmin*a+eps)),1);
        q = q + 1;
        %[q sum(nconts) a]
    end
    
    contvals = zeros(1,sum(nconts));
    
    for i = 1:nstep
    
        valsstep = fluxsteps(i) + (fluxsteps(2)-fluxsteps(1)) * ((1:nconts(i))-1)./nconts(i);
        
        indend = sum(nconts(1:i));
        indstart = indend - nconts(i)+1;
        
        contvals(indstart:indend) = valsstep;
    
    end
    

    
end
        
        