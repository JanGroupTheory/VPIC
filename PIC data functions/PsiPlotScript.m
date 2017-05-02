q = 0;

for slicenum = 1:15:91
    
    q = q + 1;
    
    load(['Psi_' num2str(slicenum) '.mat'])
    load(['Ppar_' num2str(slicenum) '.mat'])
    load(['Pperp1_' num2str(slicenum) '.mat'])
    load(['Pperp2_' num2str(slicenum) '.mat'])
    load(['ne_' num2str(slicenum) '.mat'])
    load(['magB_' num2str(slicenum) '.mat'])
    
    Cmat = contourc(Psi,[7.5 7.5]);
    
    npts = Cmat(2,1);
    
    contz = Cmat(1,2:npts+1);
    contx = Cmat(2,2:npts+1);
    
    sfac = 20;
    
    PsiVals = squeeze(smooth(interp2(Psi',contx,contz),sfac));
    Ppvals = squeeze(smooth(interp2(Ppar',contx,contz),sfac));
    Ptvals = squeeze(smooth(interp2((Pperp1'+Pperp2')/2,contx,contz),sfac));
    Bmvals = squeeze(smooth(interp2(magb',contx,contz),sfac));
    nevals = squeeze(smooth(interp2(ne',contx,contz),sfac));
    dlvals = sqrt((contz - circshift(contz,[0,1])).^2+(contx - circshift(contx,[0,1])).^2);
    
    Narr(q) = sum(nevals.*dlvals'./Bmvals);
    Larr(q) = sum(dlvals);
    Varr(q) = sum(dlvals'./Bmvals);
    
    
    %figure(100*q + 1),clf,plot(PsiVals)
%     figure(100*q + 2),clf,plot(Ppvals)
%     figure(100*q + 3),clf,plot(Ptvals)
%     figure(100*q + 4),clf,plot(Bmvals)
     %figure(5),plot(cumsum(dlvals),nevals)
    
     %hold on
end
%hold off
figure(1)