function [favg, Psivals] = loopintegral(vecarrs, barrs, varargin)

    scalfac = pi*sqrt(2)/15;

    nargs = length(varargin);

    if(nargs == 1)
        Cmat = varargin{1};
    elseif(nargs == 2)
        Cmat = contourc(varargin{1}, varargin{2});
    else
        'Invalid number of arguments.'
        favg = 0;
        Psivals = 0;
        return
    end
    
    i = 1;
    q = 1;
    
    while(i <= length(Cmat))
        
        contlen = Cmat(2,i);
       
        contx = Cmat(1,(i+1):(i+contlen));
        contz = Cmat(2,(i+1):(i+contlen)); 

        if((contx(1) ~= contx(end))||(contz(1) ~= contz(end)))
            favg(:,q) = [NaN NaN];
            Psivals(q) = Cmat(1,i);
            
        else

            vecx = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            vecy = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            vecz = interp2(vecarrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            bx = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            by = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);
            bz = interp2(barrs{1},(contz + circshift(contz,[0,1]))/2,(contx + circshift(contx,[0,1]))/2);

            dcontx = contx - circshift(contx,[0,1]);
            dcontz = contz - circshift(contz,[0,1]);

            if(isempty(find((dcontx.*bx + dcontz.*bz) > 0,1)))    
                dconty = -sqrt(dcontx.^2+dcontz.^2).*by./sqrt(bx.^2+bz.^2);
            else
                dconty = sqrt(dcontx.^2+dcontz.^2).*by./sqrt(bx.^2+bz.^2);
            end

            PhPar = cumsum(vecx.*dcontx + vecy.*dconty + vecz.*dcontz)*scalfac;

            favg(:,q) = [PhPar(end),max(PhPar) - min(PhPar)];
            Psivals(q) = Cmat(1,i);
        end
        
        q = q + 1;
        i = i + 1 + contlen;
   end
    
    
end