function [intval] = interplist(fvals,xarr,yarr,xval,yval)

   W = exp(-((xarr-xval).^2 + (yarr-yval).^2)/2);
   
   intval = (W*fvals')/sum(W);