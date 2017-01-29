function [] = makemagdata(timeslices,folderpath,vecfile,vecvar,coordnames)

    ndim = length(coordnames);
    nvec = length(vecfile);

    for i = timeslices
        
        slicestr = ['_' num2str(i) '.mat'];
        magcell = cell(1,nvec);
        
        parfor j = 1:nvec

            if(~exist([folderpath vecfile{j} 'mag' slicestr],'file'))

                vmag2 = 0;
                
                for k = 1:ndim
                    vtemp = load([folderpath vecfile{j} coordnames{k} slicestr]);
                    vmag2 = vmag2 + vtemp.([vecvar{j} coordnames{k}]).^2;
                end
                
                magcell{j} = sqrt(vmag2);
                
            end
        end
        
        for j = 1:nvec

            if(~exist([folderpath 'mag' vecfile{j} slicestr],'file'))
                eval(['mag' vecvar{j}  '= magcell{j};']);
                save([folderpath 'mag' vecfile{j} slicestr],['mag' vecvar{j}]);
            end
        end
            
    end
    
end