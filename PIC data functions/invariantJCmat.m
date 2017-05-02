function [invJ] = invariantJCmat(KE, mu, exarr, ezarr, magb, Cmat)

    if((Cmat(1,2) ~= Cmat(1,end))||(Cmat(2,2) ~= Cmat(2,end)))
        'Contour not closed.'
        return
    end

    npts = Cmat(2,1);
    
    contx = Cmat(1,2:end-1);
    conty = Cmat(2,2:end-1);
    
    dcontx = circshift(contx,[0,1]) - contx;
    dconty = circshift(conty,[0,1]) - conty;
    
    exvals = continterp(exarr, contx, conty);
    ezvals = continterp(ezarr, contx, conty);
    mbvals = continterp(magb, contx, conty);
    
    smwd = floor(npts/80);
    
    if(smwd == 0)
        Fsmooth = 1;
    else
        Fsmooth = exp(-(-2*smwd:2*smwd).^2/(2*smwd^2))/sum(exp(-(-2*smwd:2*smwd).^2/(2*smwd^2)));
    end
    
    mPhpar = conv(cumsum((circshift(contx,[0,1]) - contx).*(circshift(exvals,[0,1]) + exvals) ...
        + (circshift(conty,[0,1]) - conty).*(circshift(ezvals,[0,1]) - ezvals))/2,Fsmooth,'same');
    
    vp2 = KE - mu*mbvals + mPhpar;
    
   [minE, minind] = min(vp2);
    
    if(minE > 0)
        'Passing orbit.'
        return
    end
    
    vp2s = [vp2(minind:end) vp2(1:(minind - 1)) + mPhpar(end)-mPhpar(1)]-mPhpar(minind);
    
    invJ = 1;
    
    
    
    
end
    
    