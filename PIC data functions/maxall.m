function [maxval] = maxall(arr)

    maxval = arr;
    
    while( length(maxval) > 1)
        
        maxval = max(maxval);
        
    end
    
end