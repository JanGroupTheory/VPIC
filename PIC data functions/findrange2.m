function [range] = findrange2(expression,timesteps,varargin)

    [varnames, filenames] = findvarnames(expression);
    
    load GlobalNames.mat
    
    %tic
    if(length(varnames) == 1)
        load sliceminmax.mat
        maxarrt = maxvals.(varnames{1});
        maxarr = maxarrt(timesteps);
        minarrt = minvals.(varnames{1});
        minarr = minarrt(timesteps);
        
    else
        
        maxarr = zeros(size(timesteps));
        minarr = zeros(size(timesteps));
        q = 0;

        for slicenum = timesteps

            q = q + 1;

    %        if(mod(slicenum,10) == 9)
    %            toc
    %            slicenum + 1
    %        end

            slicestr = ['_' num2str(slicenum) '.mat'];

            for i = 1:length(varnames)
                temp = [folderpath filenames{i} slicestr];
                if(exist(temp,'file') ~= 0)
                    load(temp)
                else
                    ['Could not find ' filenames{i} slicestr '. Terminating.']
                    range = [];
                    return
                end
            end
            
            eval(['expressval =' expression ';'])
            maxarr(q) = maxall(expressval);
            minarr(q) = minall(expressval);
            
        end
        
    end
    
    if(~isempty(varargin))
        if(strcmp(varargin{1}, 'mean'))
            range(1) = mean(minarr);
            range(2) = mean(maxarr);
        else
            range(1) = min(minarr);
            range(2) = max(maxarr);
        end 
    else
        range(1) = min(minarr);
        range(2) = max(maxarr);
    end
    
    if((length(varargin) > 1) && strcmp(varargin{2},'zero'))
            range(1) = min(range(1),0);
            range(2) = max(range(2),0);
    end
end