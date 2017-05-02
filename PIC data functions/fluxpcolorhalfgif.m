function [] = fluxpcolorhalfgif(titlename, expression, timesteps, varargin)
    
    tic

    load GlobalNames.mat
    load pcolorthirdCmat.mat
    
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
    
    [~, filenames] = findvarnames(expression);
    
    nstep = length(timesteps);
    
    
    for j = 1:nstep
       
        slicenum = timesteps(j);
       
        if(mod(slicenum,10) == 9)
            toc
            slicenum + 1
        end
        
        slicestr = ['_' num2str(slicenum) '.mat'];
       
        for i = 1:length(filenames)
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
        
        [a , b] = meshgrid(-10:10,-10:10);
        smwd = 4;
        
        Fsmooth = exp(-(a.^2+b.^2)/(2*smwd^2))/sumall(exp(-(a.^2+b.^2)/(2*smwd^2)));
        
        figure(99)
        clf
        pcolor(conv2(expvals,Fsmooth,'same'))
        title([titlename sprintf('\nt = ') num2str(slicenum)])
        shading interp
        colorbar
        colormap jet
        caxis(range)
        axis([0 xmax  floor(ymax/3) ceil(ymax*2/3)])
        hold on
        plotCmat(PsiCmat.(['Cmat' num2str(timesteps(j))]) ,[0 xmax ; floor(ymax/3) ceil(ymax*2/3)]);
        hold off
        
        drawnow
        frame = getframe(99);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
       
        if j == 1;
             imwrite(A,map,[titlename ' Movie.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
        else
             imwrite(A,map,[titlename ' Movie.gif'],'gif','WriteMode','append','DelayTime',0.1);
        end
        
    end
    
    toc
    
end
    
        