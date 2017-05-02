function [invJ] = invariantJCmat3(phasespace, earrs, barrs, magb, psiarr)

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

    if((max(contmat(2,:)) == psisize(2))||(max(contmat(1,:)) == psisize(1))||(min(contmat(2,:)) == 0)||(min(contmat(1,:)) == 0))
        contclose = 0;
    else
        contclose = 1;
    end

    npts = length(contmat);
    
    contxmid = (contmat(2,1:end-1) + contmat(2,2:end))/2;
    contzmid = (contmat(1,1:end-1) + contmat(1,2:end))/2;
    
    dcontx = contmat(2,2:end) -contmat(2,1:end-1);
    dcontz = contmat(1,2:end) -contmat(1,1:end-1);

%   [dcontx, dcontz] = cleandiff(contx - circshift(contx,[0,1]),contz - circshift(contz,[0,1]),psisize(1),psisize(2));

    smwd = 10;
    [a , b] = meshgrid(-2*smwd:2*smwd,-2*smwd:2*smwd);
    Fsmooth = exp(-(a.^2+b.^2)/(2*smwd^2))/sumall(exp(-(a.^2+b.^2)/(2*smwd^2)));
%   Fsmooth = 1;

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
%    lvals = cumsum(dlvals);
    
%     smwd = 0; %floor(npts/100);
%     
%     if(smwd == 0)
%         Fsmooth = 1;
%     else
%         Fsmooth = exp(-(-2*smwd:2*smwd).^2/(2*smwd^2))/sum(exp(-(-2*smwd:2*smwd).^2/(2*smwd^2)));
%     end
    
    mPhpar = [0 cumsum((dcontx.*exvals + dconty.*eyvals + dcontz.*ezvals)*scalfac)];
       
%    edl2 = (dcontx.*(circshift(exvals,[0,1]) + exvals) + ...
%            dcontz.*(circshift(ezvals,[0,1]) + ezvals))/2*scalfac;
    
    vp2 = U2m - mu2m*mbvals + ((mPhpar(1:end-1)+mPhpar(2:end))/2 - mPhpar(minind(k)));
    
%     figure(1);plot(lvals,exvals,lvals,ezvals,lvals,edl./sqrt(dcontx.^2+dcontz.^2));
%     figure(2);pcolor(conv2(earrs{1},Fsmooth,'same'));colormap jet;shading interp;colorbar
%               hold on; plot(contx(1:3000),contz(1:3000),'k'); hold off;
%               axis([1000 1400 500 1900])
%     figure(3);pcolor(conv2(earrs{3},Fsmooth,'same'));colormap jet;shading interp;colorbar
%               hold on; plot(contx(1:3000),contz(1:3000),'k'); hold off;
%               axis([1000 1400 500 1900])
    

    
    if(contclose)
        negind = find(vp2 < 0,1);
        if(isempty(negind))
            invJ = sum(sqrt(vp2).*dlvals)*scalfac;
        else
            vp2c = circshift(vp2,[0,1-negind]);
            dlc = circshift(dlvals,[0,1-negind]);
            
            knew = mod(minind(k)-(negind-1),npts);
    
            bouncehigh = knew + find(vp2c(knew:end)<0,1) - 2;
            bouncelow = knew - find(vp2c(knew:-1:1) < 0,1) + 2;
    
            invJ = sum(sqrt(vp2c(bouncelow:bouncehigh)).*dlc(bouncelow:bouncehigh))*scalfac;
        end
    else
        
        knew = minind(k);
        
        bouncehigh = knew + find(vp2(knew:end)<0,1) - 1;
        bouncelow  = knew - find(vp2(knew:-1:1) < 0,1) + 1;
    
        if((bouncelow < 1) || isempty(bouncehigh))
            invJ = Inf;
        else
            invJ = sum(sqrt(vp2(bouncelow:bouncehigh)).*dlvals(bouncelow:bouncehigh))*scalfac;
        end
    end
    
    'done?'
    
end
