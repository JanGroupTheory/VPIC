function [vecpar, xvals, zvals] = vparvals(vecarrs, barrs, varargin)

    scalfac = pi*sqrt(2)/15;

    nargs = length(varargin);

    if(nargs == 1)
        Cmat = varargin{1};
    elseif(nargs == 2)
        Cmat = contourc(varargin{1}, varargin{2});
    else
        'Invalid number of arguments.'
        vecpar = 0;
        xvals = 0;
        zvals = 0;
        return
    end
    
    i = 1;
    q = 1;
    vecpar = [];
    xvals = [];
    zvals = [];
    
    while(i <= length(Cmat))
        
        contlen = Cmat(2,i);
       
        contx = Cmat(1,(i+1):(i+contlen));
        contz = Cmat(2,(i+1):(i+contlen)); 

        if(0) %(contx(1) ~= contx(end))||(contz(1) ~= contz(end)))
            vecpar = [vecpar NaN*contx];
            
        else

            vecx = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            vecy = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            vecz = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            bx = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            by = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            bz = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);

            dcontx = contx - circshift(contx,[0,1]);
            dcontz = contz - circshift(contz,[0,1]);
            
            ap = find((dcontx.*bx + dcontz.*bz) =< 0);
            pp = find((dcontx.*bx + dcontz.*bz) > 0);
            
            dconty = zeros(size(dcontx));
            
            dconty(ap) = -sqrt(dcontx(ap).^2+dcontz(ap).^2).*by(ap)./sqrt(bx(ap).^2+bz(ap).^2);
            dconty(pp) = sqrt(dcontx(pp).^2+dcontz(pp).^2).*by(pp)./sqrt(bx(pp).^2+bz(pp).^2);

            vdl = vecx.*dcontx + vecy.*dconty + vecz.*dcontz;
            
            vecpar = [vecpar (-cumsum(vdl)+sum(vdl)/2) *scalfac];%./sqrt(dcontx.^2+dconty.^2+dcontx.^2)];

        end
        
        xvals = [xvals contx];
        zvals = [zvals contz];
        
        q = q + 1;
        i = i + 1 + contlen;
    end
    
    
end