function [] = fluxpcolorpargif(titlename, expression, timesteps, varargin)

    folderpath = '';
    
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
    Acells = cell(nstep);
    mapcells = cell(nstep);
    Cmatcells = cell(nstep);
    
    for j = 1:nstep
        Cmatcells{j} = PsiCmat.(['Cmat' num2str(timesteps(j))]);
    end
    tic
    
    parfor j = 1:nstep
       
       slicenum = timesteps(j);
       
%        if(mod(slicenum,10) == 9)
%            toc
%            slicenum + 1
%        end
        
       slicestr = ['_' num2str(slicenum) '.mat'];
       
       expvals = evalexpress(expression,filenames,folderpath,slicestr);
       
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
       plotCmat(Cmatcells{j} ,[0 xmax ; floor(ymax/3) ceil(ymax*2/3)]);
       hold off
        
       drawnow
       frame = getframe(99);
       im = frame2im(frame);
       [Acells{j},mapcells{j}] = rgb2ind(im,256); 
       
       slicenum
       
    end
    
    for j = 1:nstep
        if j == 1;
             imwrite(Acells{j},mapcells{j},[titlename ' Movie.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
        else
             imwrite(Acells{j},mapcells{j},[titlename ' Movie.gif'],'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    
    toc
    
end


function [expvals] = evalexpress(expression,filenames,folderpath,slicestr)

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
        
end
    
        