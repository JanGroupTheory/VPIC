function [Phimat,xmat,zmat] = PhiParCmat(earrs,barrs,Cmat,scalfac)

    exarr = earrs{1};
    eyarr = earrs{2};
    ezarr = earrs{3};

    bxarr = barrs{1};
    byarr = barrs{2};
    bzarr = barrs{3};

    i = 1;
    q = 1;

    Phimat = [];
    xmat = [];
    zmat = [];
    
    dims = size(exarr);

    while(i <= length(Cmat))

        contlen = Cmat(2,i);

        contvals = Cmat(:,(i+1):(i+contlen));

        exvals = interp2(exarr,contvals(1,:),contvals(2,:));
        eyvals = interp2(eyarr,contvals(1,:),contvals(2,:));
        ezvals = interp2(ezarr,contvals(1,:),contvals(2,:));

        bxvals = interp2(bxarr,contvals(1,:),contvals(2,:));
        byvals = interp2(byarr,contvals(1,:),contvals(2,:));
        bzvals = interp2(bzarr,contvals(1,:),contvals(2,:));

        [Phitemp,xtemp,ztemp] = PhiInt({exvals,eyvals,ezvals},{bxvals,byvals,bzvals},contvals,dims(2)/2,scalfac);
        
        xmat = [xmat xtemp];
        zmat = [zmat ztemp];
        Phimat = [Phimat Phitemp];

        q = q + 1;
        i = i + 1 + contlen;

    end
        
end

function [phint,xint,zint] = PhiInt(evarrs,bvarrs,cvals,midx,scalfac)

    if((cvals(1,1) == cvals(1,end)) && (cvals(2,1) == cvals(2,end)))
        
        dlperpvals = cvals - circshift(cvals,[0,1]);
        dist2vals = dlperpvals(1,:).^2 + dlperpvals(2,:).^2;
        
        if(isempty(find((dlperpvals(1,:).*bvarrs{1} + dlperpvals(2,:).*bvarrs{3}) > 0,1)))    
            dyvals = -sqrt(dlperpvals(1,:).^2+dlperpvals(2,:).^2).*bvarrs{2}./sqrt(bvarrs{1}.^2+bvarrs{3}.^2);
        else
            dyvals = sqrt(dlperpvals(1,:).^2+dlperpvals(2,:).^2).*bvarrs{2}./sqrt(bvarrs{1}.^2+bvarrs{3}.^2);
        end
        
        edl = (dlperpvals(1,:).*(circshift(evarrs{1},[0,1]) + evarrs{1}) + ...
                dyvals.*(circshift(evarrs{2},[0,1]) + evarrs{2}) + ... 
                dlperpvals(2,:).*(circshift(evarrs{3},[0,1]) + evarrs{3}))/2*scalfac;
            
        [~,cutind] = max(abs(cvals(1,:) - midx));
        
        Phpar = cumsum(circshift(edl,[0,-cutind]));
        
        xtemp = circshift(cvals(1,:),[0,-cutind]);
        ztemp = circshift(cvals(2,:),[0,-cutind]);
        phtemp = Phpar - (Phpar(end) + Phpar(1))/2;
        
        nzinds = find(circshift(dist2vals,[0,-cutind]) ~= 0);
        
        phint = phtemp(nzinds);
        xint = xtemp(nzinds);
        zint = ztemp(nzinds);
        
        
    else
        
        dlperpvals = [cvals(:,2:end) - cvals(:,1:end-1) [0;0]];
        dist2vals = dlperpvals(1,:).^2 + dlperpvals(2,:).^2;
        
        if(isempty(find((dlperpvals(1,:).*bvarrs{1} + dlperpvals(2,:).*bvarrs{3}) > 0,1)))    
            dyvals = -sqrt(dlperpvals(1,:).^2+dlperpvals(2,:).^2).*bvarrs{2}./sqrt(bvarrs{1}.^2+bvarrs{3}.^2);
        else
            dyvals = sqrt(dlperpvals(1,:).^2+dlperpvals(2,:).^2).*bvarrs{2}./sqrt(bvarrs{1}.^2+bvarrs{3}.^2);
        end
        
        edl = (dlperpvals(1,:).*(circshift(evarrs{1},[0,1]) + evarrs{1}) + ...
                dyvals.*(circshift(evarrs{2},[0,1]) + evarrs{2}) + ... 
                dlperpvals(2,:).*(circshift(evarrs{3},[0,1]) + evarrs{3}))/2*scalfac;
        
        Phpar = cumsum(edl);
        
        phtemp = Phpar - (Phpar(end) + Phpar(1))/2;
        nzinds = find(dist2vals ~= 0);
        
        phint = phtemp(nzinds);
        xint = cvals(1,nzinds);
        zint = cvals(2,nzinds);
    end

end