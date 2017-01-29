function [] = plotCmat(Cmat,varargin)

    if(~isempty(varargin))
        lims = varargin{1};
    end

    i = 1;
    
    hold on
    while(i <= length(Cmat))
        
        contlen = Cmat(2,i);
        
        contvals = Cmat(:,(i+1):(i+contlen));
        
        plot(contvals(1,:),contvals(2,:),'k');
        xlim(lims(1,:))
        ylim(lims(2,:))
        
        i = i + 1 + contlen;
      
        
    end
    hold off