function [invJ] = invariantJCmat4(phasespace, earrs, barrs, psiarr)

%    scalfac = pi*sqrt(2)/15;
    Fsmooth = 1;
%     smwd = 10;
%     [a , b] = meshgrid(-2*smwd:2*smwd,-2*smwd:2*smwd);
%     Fsmooth = exp(-(a.^2+b.^2)/(2*smwd^2))/sumall(exp(-(a.^2+b.^2)/(2*smwd^2)));
    scalfac = 1/200;

    x = phasespace(1);
    z = phasespace(2);
    vx = phasespace(3);
    vy = phasespace(4);
    vz = phasespace(5);
    
    magb = sqrt(barrs{1}.^2 + barrs{2}.^2 + barrs{3}.^2);
    
    psisize = size(psiarr);
    
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

    if((max(contmat(2,:)) == psisize(2))||(max(contmat(1,:)) == psisize(1))...
            ||(min(contmat(2,:)) == 1)||(min(contmat(1,:)) == 1))
        contclose = 0;
        'open'
    else
        'closed'
        contclose = 1;
    end
    
    if(contclose)
        contx = circshift(contmat(2,1:end),[0,-minind(k)]);
        contz = circshift(contmat(1,1:end),[0,-minind(k)]);
    
        dcontx = contx - circshift(contx,[0,1]);
        dcontz = contz - circshift(contz,[0,1]);
        
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
        
        edl = (dcontx.*(circshift(exvals,[0,1]) + exvals) + ...
                dconty.*(circshift(eyvals,[0,1]) + eyvals) + ... 
                dcontz.*(circshift(ezvals,[0,1]) + ezvals))/2*scalfac;
       
        mPhparCCW = conv(cumsum(edl),1,'same');
        mPhparCW = conv(cumsum(-edl(end:-1:1)),1,'same');
    
        vp2CCW = U2m - mu2m*mbvals + mPhparCCW;
        vp2CW = U2m - mu2m*mbvals(end:-1:1) + mPhparCW;
        vp2CCW(1) = max(0,vp2CCW(1));
        vp2CW(1) = max(0,vp2CCW(1));
        
        bounceCCW = find(vp2CCW < 0,1);
        bounceCW = find(vp2CW < 0,1);
        
        if(isempty(bounceCCW)||isempty(bounceCW))
            invJ = Inf;
            return
        end
        
        invJ = (sum(dlvals(1:bounceCCW-1).*sqrt(vp2CCW(1:bounceCCW-1))) ... 
                    +sum(dlvals(end:-1:end-bounceCW+2).*sqrt(vp2CW(1:bounceCW-1))))*scalfac;
        
    else
        npts = length(contmat);
    
        contxmid = (contmat(2,1:end-1) + contmat(2,2:end))/2;
        contzmid = (contmat(1,1:end-1) + contmat(1,2:end))/2;
    
        dcontx = contmat(2,2:end) -contmat(2,1:end-1);
        dcontz = contmat(1,2:end) -contmat(1,1:end-1);
        
        exvals = interp2(conv2(earrs{1},Fsmooth,'same'), contzmid, contxmid);
        eyvals = interp2(conv2(earrs{2},Fsmooth,'same'), contzmid, contxmid);
        ezvals = interp2(conv2(earrs{3},Fsmooth,'same'), contzmid, contxmid);
        bxvals = interp2(conv2(barrs{1},Fsmooth,'same'), contzmid, contxmid);
        byvals = interp2(conv2(barrs{2},Fsmooth,'same'), contzmid, contxmid);
        bzvals = interp2(conv2(barrs{3},Fsmooth,'same'), contzmid, contxmid);
        mbvals = interp2(conv2(magb,Fsmooth,'same'), contzmid, contxmid);

        if(isempty(find((dcontx.*bxvals + dcontz.*bzvals) > 0,1)))    
            dconty = -sqrt(dcontx.^2+dcontz.^2).*byvals./sqrt(mbvals.^2-byvals.^2);
        else
            dconty = sqrt(dcontx.^2+dcontz.^2).*byvals./sqrt(mbvals.^2-byvals.^2);
        end

        dlvals = sqrt(dcontx.^2+dconty.^2+dcontz.^2);
        
        mPhpar = [0 cumsum((dcontx.*exvals + dconty.*eyvals + dcontz.*ezvals)*scalfac)];

        vp2 = U2m - mu2m*mbvals + ((mPhpar(1:end-1)+mPhpar(2:end))/2 - mPhpar(minind(k)));
        
        knew = minind(k);
        
        if((knew == 1)||(knew > length(vp2)))
            invJ = Inf;
            return
        end
        
        bouncehigh = knew + find(vp2(knew:end)<0,1) - 2;
        bouncelow  = knew - find(vp2(knew:-1:1) < 0,1) + 2;
    
        if(isempty(bouncehigh)||isempty(bouncelow)||(bouncelow < 1))
            invJ = Inf;
        else
            invJ = sum(sqrt(vp2(bouncelow:bouncehigh)).*dlvals(bouncelow:bouncehigh))*scalfac;
        end
        
    end