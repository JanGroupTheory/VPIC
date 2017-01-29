function [] = fluxpcolorgif(titlename, expression, timesteps, varargin)

    %Creates a density plot gif of an arbitrary function of the plasma
    %variables with overlaid flux surfaces at the specified timesteps. The
    %flux surfaces are drawn from the precalculated structure created by
    %makefluxCmat.
    
    %The argument 'titlename' is a string that specifies the title of the
    %plot. This will also be used for the name of the resulting file.
    
    %The argument 'expression' is a string that specifies what should be
    %plotted. It should be given in the form of a MATLAB command that
    %includes only the variable names specified in setup.
    
    %The argument 'timesteps' is a list of the timesteps to include in the
    %gif. This should correspond to the numerical names of the timesteps in
    %the variable file names.
   
    %The optional argument specifies the range to use for the density plot.
    %If unspecified, the function will calculate the range on its own, which 
    %takes more time. If specified as 'zero', the function will ensure that
    %the calculated range includes zero. 
    
    load GlobalNames.mat 
    load pcolorCmat.mat 
    
    if(isempty(varargin)||strcmp(varargin{1},''))
        range = findrange(expression, timesteps, 'mean');
    elseif(strcmp(varargin{1}, 'zero'))
        range = findrange(expression, timesteps, 'mean','zero');
    else
        range = varargin{1};
    end
    %Finds the range to use for the plot.
    
    [varnames, filenames] = findvarnames(expression);
    %Finds all variables used in expression and their corresponding file
    %names.
    
    %tic
    
    for slicenum = timesteps
       
        %if(mod(slicenum,10) == 9)
        %    toc
        %    slicenum + 1
        %end 
            
        slicestr = ['_' num2str(slicenum) '.mat'];
        
        for i = 1:length(varnames)
            temp = [folderpath filenames{i} slicestr];
            if(exist(temp,'file') ~= 0)
                load(temp)
            else
                ['Could not find' filenames{i} '_' slicestr '. Terminating.']
                return
            end
        end
       
        eval(['pcolor(' expression ')'])
        
        clf
        title([titlename sprintf('\nt = ') num2str(slicenum)])
        shading interp
        colorbar
        colormap jet
        caxis(range)
        hold on
        plotCmat(PsiCmat.(['Cmat' num2str(slicenum)]) ,[0 2560 ; 0 1280]);
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
    
    %toc
    
end
        