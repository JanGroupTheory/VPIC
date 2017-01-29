function [favg, Cmat] = fluxsurfavg(farr1, varargin)

    nargs = length(varargin);

    if(nargs == 1)
        Cmat = varargin{1};
    elseif(nargs == 2)
        [Cmat , ~] = contour(varargin{1}, varargin{2});
    else
        'Invalid number of arguments.'
        favg = 0;
        Cmat = 0;
        return
    end
    
    i = 1;
    q = 1;
    
    if(length(farr1) == 1)
        farr = farr1*ones([ceil(max(abs(Cmat(2,:)))), ceil(max(abs(Cmat(1,:))))]);
    else
        farr = farr1;
    end
    
    while(i <= length(Cmat))
        
        contlen = Cmat(2,i);
        
        contvals = Cmat(:,(i+1):(i+contlen));
        
%         if(prod(contvals(:,1) == contvals(:,end)) ~= 1)
%             
%              'Contour not closed.'
%              fint = 0;
%              i = i + 1 + contlen;
%              continue
%         end
        
        xp = ceil(contvals(1,:));
        xm = floor(contvals(1,:));
        dx = contvals(1,:) - xm;
        yp = ceil(contvals(2,:));
        ym = floor(contvals(2,:));
        dy = contvals(2,:) - ym;
           
        fmm = farr(sub2ind(size(farr),ym,xm));
        fpm = farr(sub2ind(size(farr),yp,xm));
        fmp = farr(sub2ind(size(farr),ym,xp));
        %fpp = farr(sub2ind(size(farr),yp,xp));
        
        fvals = fmm.*(1 - dx).*(1 - dy) + fmp.*dx.*(1-dy) + fpm.*dy.*(1-dx);% + (1-dx).*(1-dy).*fpp
        
        
        if(prod(contvals(:,1) == contvals(:,end)) == 1)
            dlvals  = sqrt(sum((contvals - circshift(contvals,1,2)).^2));
            fint = sum(dlvals.*(fvals + circshift(fvals,1,2))/2);
        else
            dlvals = sqrt(sum((contvals(:,2:end) - contvals(:,1:end-1)).^2));
            fint = sum(dlvals.*(fvals(2:end) + fvals(1:end-1))/2);
        end
        
        if(sum(dlvals) == 0)
            favg(q) = meanall(fvals);
        else
            favg(q) = fint/sum(dlvals);
        end
        
        if(isnan(favg(q)))
            'God why'
        end
        
        q = q + 1;
        i = i + 1 + contlen;
   end
    
    
end