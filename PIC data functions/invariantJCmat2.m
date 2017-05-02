function [invJ] = invariantJCmat2(phasespace, earrs, barrs, magb, psiarr)

    psisize = size(psiarr);

    scalfac = pi*sqrt(2)/15;
%    scalfac = 1/200;

    x = phasespace(1);
    z = phasespace(2);
    vx = phasespace(3);
    vy = phasespace(4);
    vz = phasespace(5);
    
    psval = interp2(psiarr,z,x);
    bxval = interp2(barrs{1},z,x);
    byval = interp2(barrs{2},z,x);
    bzval = interp2(barrs{3},z,x);
    mbval = interp2(magb,z,x);
    
    Cmat = contourc(psiarr,[psval psval]);
    
    mu2m = sum(cross([vx,vy,vz],[bxval,byval,bzval]).*cross([vx,vy,vz],[bxval,byval,bzval]))/(mbval^3);
    U2m = (vx^2+vy^2+vz^2);
    
    q = 1;
    i = 1;
    
    while(i <= length(Cmat))
        contlen = Cmat(2,i);
        
        contvals{q} = Cmat(:,(i+1):(i+contlen));
        
        q = q + 1;
        i = i + 1 + contlen;
    end

    for j = 1:length(contvals)
       contmat = contvals{j};
       [mindist(j),minind(j)] = min((contmat(2,:)-x).^2 + (contmat(1,:)-z).^2);
    end
    
    [~,k] = min(mindist);
    
    contmat = contvals{k};
    
%     if((contmat(1,1) ~= contmat(1,end))||(contmat(2,1) ~= contmat(2,end)))
%         'Contour not closed.'
%         invJ = Inf;
%         return
%     end

%    npts = length(contmat);
    
    contx = circshift(contmat(2,1:end),[0,-minind(k)]);
    contz = circshift(contmat(1,1:end),[0,-minind(k)]);
    
    dcontx = contx - circshift(contx,[0,1]);
    dcontz = contz - circshift(contz,[0,1]);

%    [dcontx, dcontz] = cleandiff(contx - circshift(contx,[0,1]),contz - circshift(contz,[0,1]),psisize(1),psisize(2));

    smwd = 10;
    [a , b] = meshgrid(-2*smwd:2*smwd,-2*smwd:2*smwd);
    Fsmooth = exp(-(a.^2+b.^2)/(2*smwd^2))/sumall(exp(-(a.^2+b.^2)/(2*smwd^2)));
%    Fsmooth = 1;

%     Fsmooth =  [13    21    31    40    48    50    48    40    31    21    13
%                 21    34    50    67    79    83    79    67    50    34    21
%                 31    50    74    98   116   123   116    98    74    50    31
%                 40    67    98   130   153   162   153   130    98    67    40
%                 48    79   116   153   181   192   181   153   116    79    48
%                 50    83   123   162   192   208   192   162   123    83    50
%                 48    79   116   153   181   192   181   153   116    79    48
%                 40    67    98   130   153   162   153   130    98    67    40
%                 31    50    74    98   116   123   116    98    74    50    31
%                 21    34    50    67    79    83    79    67    50    34    21
%                 13    21    31    40    48    50    48    40    31    21    13]/10000;
    
    exvals = interp2(conv2(earrs{1},Fsmooth,'same'), contz, contx);
    eyvals = interp2(conv2(earrs{2},Fsmooth,'same'), contz, contx);
    ezvals = interp2(conv2(earrs{3},Fsmooth,'same'), contz, contx);
    bxvals = interp2(conv2(barrs{1},Fsmooth,'same'), contz, contx);
    byvals = interp2(conv2(barrs{2},Fsmooth,'same'), contz, contx);
    bzvals = interp2(conv2(barrs{3},Fsmooth,'same'), contz, contx);
    mbvals = interp2(conv2(magb,Fsmooth,'same'), contz, contx);
    
    if(isempty(find((dcontx.*bxvals + dcontz.*bzvals) > 0,1)))    
        dconty = -sqrt(dcontx.^2+dcontz.^2).*byvals./sqrt(mbvals.^2-byvals.^2);
    else
        dconty = sqrt(dcontx.^2+dcontz.^2).*byvals./sqrt(mbvals.^2-byvals.^2);
    end
    
    dlvals = sqrt(dcontx.^2+dconty.^2+dcontz.^2);
%    lvals = cumsum(dlvals);
    
%     smwd = 0; %floor(npts/100);
%     
%     if(smwd == 0)
%         Fsmooth = 1;
%     else
%         Fsmooth = exp(-(-2*smwd:2*smwd).^2/(2*smwd^2))/sum(exp(-(-2*smwd:2*smwd).^2/(2*smwd^2)));
%     end
    
    edl = (dcontx.*(circshift(exvals,[0,1]) + exvals) + ...
           dconty.*(circshift(eyvals,[0,1]) + eyvals) + dcontz.*(circshift(ezvals,[0,1]) + ezvals))/2*scalfac;
       
%    edl2 = (dcontx.*(circshift(exvals,[0,1]) + exvals) + ...
%            dcontz.*(circshift(ezvals,[0,1]) + ezvals))/2*scalfac;
    
    mPhparCCW = conv(cumsum(edl),1,'same');
    mPhparCW = conv(cumsum(-edl(end:-1:1)),1,'same');
    
    vp2CCW = U2m - mu2m*mbvals + mPhparCCW;
    vp2CW = U2m - mu2m*mbvals(end:-1:1) + mPhparCW;
    vp2CCW(1) = max(0,vp2CCW(1));
    vp2CW(1) = max(0,vp2CCW(1));
    
%     figure(1);plot(lvals,exvals,lvals,ezvals,lvals,edl./sqrt(dcontx.^2+dcontz.^2));
%     figure(2);pcolor(conv2(earrs{1},Fsmooth,'same'));colormap jet;shading interp;colorbar
%               hold on; plot(contx(1:3000),contz(1:3000),'k'); hold off;
%               axis([1000 1400 500 1900])
%     figure(3);pcolor(conv2(earrs{3},Fsmooth,'same'));colormap jet;shading interp;colorbar
%               hold on; plot(contx(1:3000),contz(1:3000),'k'); hold off;
%               axis([1000 1400 500 1900])
    
    bounceCCW = find(vp2CCW < 0,1);
    bounceCW = find(vp2CW < 0,1);
    
    if(isempty(bounceCCW)||isempty(bounceCW))
%        'Passing orbit.'
        invJ = Inf;
        return
    end
    
    invJ = (sum(dlvals(1:bounceCCW-1).*sqrt(vp2CCW(1:bounceCCW-1))) ... 
        +sum(dlvals(end:-1:end-bounceCW+2).*sqrt(vp2CW(1:bounceCW-1))))*scalfac;
  
%    'done?'
    
end

function [cdx, cdz] = cleandiff(dx,dz,xsize,zsize)
    
    cdx = dx;
    cdz = dz;
    
    if(dx(1)^2 > 2)
        if((dx(1) - xsize)^2 < 2)
            cdx(1) = dx(1) - xsize;
        elseif((dx(1) + 2400)^2 < 2)
            cdx(1) = dx(1) + xsize;
        else
            'Contour not closed.'
            cdx = Inf;
            return
        end
    end
    
    if(dx(end)^2 > 2)
        if((dx(end) - xsize)^2 < 2)
            cdx(end) = dx(end) - xsize;
        elseif((dx(1) + xsize)^2 < 2)
            cdx(end) = dx(end) + xsize;
        else
            'Contour not closed.'
            cdx = Inf;
            return
        end
    end
    
    if(dz(end)^2 > 2)
        if((dz(end) - zsize)^2 < 2)
            cdz(end) = dz(end) - zsize;
        elseif((dx(1) + zsize)^2 < 2)
            cdz(end) = dz(end) + zsize;
        else
            'Contour not closed.'
            cdx = Inf;
            return
        end
     end
    
    if(dz(1)^2 > 2)
        if((dz(1) - zsize)^2 < 2)
            cdz(1) = dz(1) - zsize;
        elseif((dx(1) + zsize)^2 < 2)
            cdz(1) = dz(1) + zsize;
        else
            'Contour not closed.'
            cdx = Inf;
            return
        end
    end

end