function [meanval] = meanall(arr)

    meanval = arr;
    
    while( length(meanval) > 1)
        
        meanval = mean(meanval);
        
    end
    
end