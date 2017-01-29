function [col] = colormapinterp( value, range, colmap)

    cmap = colormap(colmap);
    nc = size(cmap);
    
    rescal = (value - range(1))/(range(2)-range(1))*(nc(1)-1) + 1;
    
    if(rescal >= nc(1))
        col =  cmap(end,:);
        return
    elseif(rescal <= 1)
        col = cmap(1,:);
        return
    end
    
    ind = floor(rescal);
    fracval = rescal - ind;
    
    if((ind < 1)||isnan(ind))
        'God dammit.'
    end
    
    col = cmap(ind,:)*fracval + cmap(ind+1,:)*(1-fracval);
    
end
    
        
    