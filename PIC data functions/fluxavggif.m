function [] = fluxavggif(titlename, expression, timesteps, varargin)

    scalfilenames = {'Psi','Ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2'};
    scalvarnames = {'Psi','ne','Ppar','Pparp1','Pparp2','Pp1p2','Pperp1','Pperp2'};
    vecfilenames = {'B','J','E','Ue'};
    vecvarnames = {'b','j','e','ue'};
    symtensnames = {'P'};
    coordnames = {'x','y','z'};
    ndim = length(coordnames);
    
    folderpath = 'D:\IslandCoalescenceTimeSlices\';
    
    tic
    
    for slicenum = timesteps
        
        if(mod(slicenum,2) == 1)
           toc
           slicenum + 1
        end
        
        slicestr = ['_' num2str(4*slicenum + 1) '.mat'];
        
        temp = [folderpath 'Psi' slicestr];
        if(exist(temp,'file') ~=0)
            load(temp)
        else
            ['Could not find Psi_' slicestr '. Terminating.']
            return
        end
        
        if(slicenum == timesteps(1))
            if(isempty(varargin))
                contvals = findcontspacing(Psi);
            else
                contvals = varargin{1};
            end
            frange = findavgrange(expression, timesteps, contvals);
        end
        
        
        
        for i = 1:length(scalvarnames)
            if(~isempty(strfind(expression,scalvarnames{i})))
                temp = char(strcat(folderpath,scalfilenames{i},slicestr));
                if(exist(temp,'file') ~= 0)
                    load(temp)
                else
                    ['Could not find ' scalfilenames{i} '_' slicestr '. Terminating.']
                    return
                end
            end  
        end

        for i = 1:length(vecvarnames)
            for coord = coordnames
                vecname = strcat(vecvarnames{i},coord);
                filename = strcat(vecfilenames{i},coord);
                if(~isempty(strfind(expression,vecname)))
                    temp = char(strcat(folderpath,filename,slicestr));
                    if(exist(temp,'file') ~= 0)
                        load(temp)
                    else
                        ['Could not find' filename '_' slicestr '. Terminating.']
                        return
                    end
                end  
            end
        end

        eval(['contavgplot(' expression ', Psi'', contvals, frange)']);
        title([titlename sprintf('\nt = ') num2str(4*slicenum + 1)])
        
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); 
       
        if slicenum == timesteps(1);
             imwrite(A,map,[titlename ' Movie.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
        else
             imwrite(A,map,[titlename ' Movie.gif'],'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    
    toc
    
end