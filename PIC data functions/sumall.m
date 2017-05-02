function [sumval] = sumall(arr)

    sumval = arr;
    
    while( length(sumval) > 1)
        
        sumval = sum(sumval);
        
    end
    
end