function [] = fluxpcolorgif(titlename, expression, timesteps, varargin)
    
    load GlobalNames.mat
    load pcolorCmat.mat
    
    if(isempty(varargin)||strcmp(varargin{1},''))
        range = findrange2(expression, timesteps, 'mean');
    elseif(strcmp(varargin{1}, 'zero'))
        range = findrange2(expression, timesteps, 'mean','zero');
    else
        range = varargin{1};
    end
    
    if(isempty(range))
        ['Invalid range for ' expression '. Terminating.']
        return
    end
    
    [varnames, filenames] = findvarnames(expression);
    
    tic
    
    for slicenum = timesteps
        
       if(mod(slicenum,10) == 9)
           toc
           slicenum + 1
       end
        
        slicestr = ['_' num2str(slicenum) '.mat'];
        
        for i = 1:length(varnames)
            temp = [folderpath filenames{i} slicestr];
            if(exist(temp,'file') ~= 0)
                load(temp)
            else
                ['Could not find ' filenames{i} slicestr '. Terminating.']
                return
            end
        end
        
        expvals = eval(expression);
        
        [ymax,xmax] = size(expvals);
        
        clf
        pcolor(expvals)
        title([titlename sprintf('\nt = ') num2str(slicenum)])
        shading interp
        colorbar
        colormap jet
        caxis(range)
        hold on
        plotCmat(PsiCmat.(['Cmat' num2str(slicenum)]) ,[0 xmax ; 0 ymax]);
        hold off
        
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
        