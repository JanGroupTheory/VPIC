function [] = maketempdata(timeslices,folderpath)

    for i = timeslices
        
        slicestr = ['_' num2str(i) '.mat'];
        
        if(~exist([folderpath 'Tperp1' slicestr],'file'))
            if(~exist('ns','var'))
                ns = load([folderpath 'Ne' slicestr]);
            end
            Pts1 = load([folderpath 'Pperp1' slicestr]);
            Tperp1 = Pts1.Pperp1./ns.ne;
            save([folderpath 'Tperp1' slicestr],'Tperp1');
        end
        
        if(~exist([folderpath 'Tperp1' slicestr],'file'))
            if(~exist('ns','var'))
                ns = load([folderpath 'Ne' slicestr]);
            end
            Pts2 = load([folderpath 'Pperp2' slicestr]);
            Tperp2 = Pts2.Pperp2./ns.ne;
            save([folderpath 'Tperp2' slicestr],'Tperp2');
        end
        
        if(~exist([folderpath 'Tperp1' slicestr],'file'))
            if(~exist('ns','var'))
                ns = load([folderpath 'Ne' slicestr]);
            end
            Pps = load([folderpath 'Ppar' slicestr]);
            Tpar = Pps.Ppar./ns.ne;
            save([folderpath 'Tpar' slicestr],'Tpar');
        end
        

    end
    
end