function [] = continterpplot(funcvals, fluxvals, contvals, varargin)

    tic
    
    if(isempty(contvals))
        [Cmat, ~] = contour(fluxvals, findcontspacing(fluxvals));
    elseif(min(size(contvals)) == 1)
        [Cmat, ~] = contour(fluxvals, contvals);
    elseif(min(size(contvals)) == 2)
        Cmat = contvals;
    else
        'Invalid value for contvals'
        return
    end
    
    
    favgvals = fluxsurfavg(funcvals,Cmat);
    
    if(isempty(varargin))
        frange = [minall(favgvals),maxall(favgvals)];
    else
        frange = varargin{1};
    end
    
    i = 1;
    q = 1;
    
    contfluxvals = zeros(size(favgvals));
    contptsarr = cell(size(favgvals));
    
    while(i <= length(Cmat))
        
        contfluxvals(q) = Cmat(1,i);
        contlen = Cmat(2,i);
        
        contptsarr{q} = Cmat(:,(i+1):(i+contlen));
        
        i = i + contlen + 1;
        q = q + 1;
    end
    
%     downsamp = 8;
%     
%     funcavg = zeros(floor(size(funcvals)/downsamp));

    funcavg = zeros(size(funcvals));

    toc
    
    for i = 1:numel(funcavg)
        
        psival = fluxvals(i);
        [pti,ptj] = ind2sub(size(funcavg),i);
        
        psilow = max(contfluxvals(contfluxvals < psival));
        psihigh = min(contfluxvals(contfluxvals > psival));
        if(isempty(psilow))
            psilow = psihigh;
        end
        if(isempty(psihigh))
            psihigh = psilow;
        end
       
        
        qlow = find(contfluxvals == psilow);
        qhigh = find(contfluxvals == psihigh);
        
        if(length(qlow) > 1)
            
            d2min = zeros(size(qlow));
            for j = 1:length(qlow)
                contpts = contptsarr{qlow(j)};
                dxvals = contpts(1,:) - pti;
                dyvals = contpts(2,:) - ptj;
                d2min(j) = min(sum(dxvals.^2 + dyvals.^2));
            end
            
            [~, js] = min(d2min);
            
            flow = favgvals(qlow(js));
        else
            flow = favgvals(qlow);
        end
        if(length(qhigh) > 1)
            
            d2min = zeros(size(qhigh));
            for j = 1:length(qhigh)
                contpts = contptsarr{qhigh(j)};
                dxvals = contpts(1,:) - pti;
                dyvals = contpts(2,:) - ptj;
                d2min(j) = min(dxvals.^2 + dyvals.^2);
            end
            
            [~, js] = min(d2min);
            
            fhigh = favgvals(qhigh(js));
        else
            fhigh = favgvals(qhigh);
        end        
        
        if(psilow == psihigh)
            funcavg(i) = fhigh;
        else
            funcavg(i) = (fhigh * (psival - psilow) + flow * (psihigh - psival))/(psihigh - psilow);
        end
        
%         wtvals = exp(-1.*sum((ptsarray(1:2,:) - repmat([pti ; ptj],[1,length(ptsarray)])).^2));
%         
%         funcavg(i) = sum( ptsarray(3,:) .* wtvals)/sum(wtvals);
        
%         if((mod(ptj,160) == 0) &&(pti == 1280))
%             [num2str(pti) ' ' num2str(ptj) ' ' num2str(round(100 * i/numel(funcavg),2)) '%']
%             toc
%         end
        
    end
    
    toc
    
    pcolor(funcavg)
    shading interp
    colorbar
    colormap jet
    caxis(frange)
    
end