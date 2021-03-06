function [varnames, filenames] = findvarnames(express)

    load GlobalNames.mat
    
    varnames = cell(0);
    filenames = cell(0);
    
    for j = 1:length(allvarnames)
        if(regexp(express,['\<' allvarnames{j} '\>']))
            varnames = cat(2,varnames,allvarnames(j));
            filenames = cat(2,filenames,allfilenames(j));
        end
    end
        
        

