function fvals = continterp(farr, contx, conty)

    xp = ceil(contx);
    xm = floor(contx);
    dx = contx - xm;
    yp = ceil(conty);
    ym = floor(conty);
    dy = conty - ym;

    fmm = farr(sub2ind(size(farr),ym,xm));
    fpm = farr(sub2ind(size(farr),yp,xm));
    fmp = farr(sub2ind(size(farr),ym,xp));
    fpp = farr(sub2ind(size(farr),yp,xp));

    fvals = fmm.*(1 - dx).*(1 - dy) + fmp.*dx.*(1-dy) + fpm.*dy.*(1-dx) + (1-dx).*(1-dy).*fpp;
    
end