function [] = fluxinterpgif(titlename, expression, timesteps, varargin)
    
    folderpath = 'D:\IslandCoalescenceTimeSlices\';
    
    tic
    
    frange = findavgrange2(expression, timesteps);
    [varnames, filenames] = findvarnames(expression);
    
    load FluxInterpCmat.mat
    
    
    for slicenum = timesteps
        
        if(mod(slicenum, 10) == 1)
           toc
           slicenum
        end
        
        Cmat = PsiCmat.(['Cmat' num2str(slicenum)]);
        
        slicestr = ['_' num2str(4*slicenum + 1) '.mat'];
        
        temp = [folderpath 'Psi' slicestr];
        if(exist(temp,'file') ~=0)
            load(temp)
        else
            ['Could not find Psi_' slicestr '. Terminating.']
            return
        end

        for i = 1:length(varnames)
            temp = [folderpath filenames{i} slicestr];
            if(exist(temp,'file') ~= 0)
                load(temp)
            else
                ['Could not find' filenames{i} '_' slicestr '. Terminating.']
                return
            end
        end

        eval(['continterpplot(' expression ', Psi, Cmat, frange)']);
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