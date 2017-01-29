function [minval] = minall(arr)

    minval = arr;
    
    while( length(minval) > 1)
        
        minval = min(minval);
        
    end
    
end