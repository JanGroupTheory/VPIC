    scalnames = {'Psi','Ne'};
    vecnames = {'B','J','E','Ue'};
    symtensnames = {'P'};
    coordnames = {'x','y','z'};
    ndim = length(coordnames);
    othernames = {'Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2'};

    slicenum = input('Which slice to read? ');

    slicestr = ['_' num2str(slicenum) '.mat'];
    
    folderpath = 'D:\ForceFreeIslands\';

    for sname = scalnames
        temp = char(strcat(folderpath,sname,slicestr));
            if(exist(temp,'file') ~= 0)
               load(temp)
               %['Loaded ' char(strcat(sname,'.'))]
            else
               ['Could not find ' char(strcat(sname,'.'))] 
            end
    end
    
    for vname = vecnames
        for i = 1:ndim
            temp = char(strcat(folderpath,vname,coordnames(i),slicestr));
            if(exist(temp,'file') ~= 0)
                load(temp)
                %['Loaded ' char(strcat(vname,coordnames(i),'.'))]
            else
                ['Could not find ' char(strcat(vname,coordnames(i),'.'))]
            end
        end
    end
    
    for sname = symtensnames
        for i = 1:ndim
            for j = i:ndim
                temp = char(strcat(folderpath,sname,coordnames(i),coordnames(j),slicestr));
                if(exist(temp,'file') ~= 0)
                    load(temp)
                    %['Loaded ' char(strcat(sname,coordnames(i),coordnames(j),'.'))]
                else
                    ['Could not find ' char(strcat(sname,coordnames(i),coordnames(j),'.'))]
                end
            end
        end
    end
    
    for oname = othernames
        temp = char(strcat(folderpath,oname,slicestr));
        if(exist(temp,'file') ~= 0)
            load(temp)
            %['Loaded ' char(strcat(oname,'.'))]
        else
            ['Could not find ' char(strcat(oname,'.'))]
        end
    end
    
    
    
%    load(['D:\Island_Coalescence_Time_Slices\' 'Bx' slicestr])
%     load ['By_' slicestr '.mat']
%     load ['Bz_' slicestr '.mat']
%     load 'Ex_' slicestr '.mat'
%     load 'Ey_' slicestr '.mat'
%     load 'Ez_' slicestr '.mat'
%     load 'Uex_' slicestr '.mat'
%     load 'Uey_' slicestr '.mat'
%     load 'Uez_' slicestr '.mat'
%     load 'Jx_' slicestr '.mat'
%     load 'Jy_' slicestr '.mat'
%     load 'Jz_' slicestr '.mat'
%     load 'Ne_' slicestr '.mat'
%     load 'Psi_' slicestr '.mat'
%     load 'Pxx_' slicestr '.mat'
%     load 'Pxy_' slicestr '.mat'
%     load 'Pxz_' slicestr '.mat'
%     load 'Pyy_' slicestr '.mat'
%     load 'Pyz_' slicestr '.mat'
%     load 'Pzz_' slicestr '.mat'
%     load 'Pperp1_' slicestr '.mat'
%     load 'Pperp2_' slicestr '.mat'
%     load 'Pp1p2_' slicestr '.mat'
%     load 'Pparp1_' slicestr '.mat'
%     load 'Pparp2_' slicestr '.mat'
%     load 'Ppar_' slicestr '.mat'
%     
%    