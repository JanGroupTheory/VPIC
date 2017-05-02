function [] = plot2list(farr,xarr,zarr,varargin)

    if(isempty(varargin))
        frange = [min(farr) max(farr)];
    else
        frange = varargin{1}; 
    end
    
    clf
    hold on
    pointinds = farr <= frange(1);
    col = colormapinterp(frange(1),frange,'jet');
    plot(xarr(pointinds),zarr(pointinds),'.','Color',col)
    
    for i = 1:100
        
        fslicemax = (frange(2) - frange(1))* i/100 + frange(1);
        fslicemin = (frange(2) - frange(1))* (i-1)/100 + frange(1);
        
        pointinds = (farr <= fslicemax) & (farr > fslicemin);
        
        col = colormapinterp((fslicemin + fslicemax)/2,frange,'jet');
        
        plot(xarr(pointinds),zarr(pointinds),'.','Color',col);
        
    end
    
    pointinds = farr > frange(2);
    col = colormapinterp(frange(2),frange,'jet');
    plot(xarr(pointinds),zarr(pointinds),'.','Color',col)
    
    hold off
    
    axis([min(xarr) max(xarr) min(zarr) max(zarr)])
end