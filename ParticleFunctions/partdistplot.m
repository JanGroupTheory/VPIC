function partdistplot(directory, domains, zs, varargin)
cd(directory)
load('partsetup.mat')
cntx=0; 
% zs = 19.5:.125:25;
Vlim=1.5;
colbnds=[-2.5 -.35];
colbndsperp = colbnds;

if nargin>3
    partype = varargin{1};
else
    partype = 'e';
end

Xps = zeros(length(domains),1); 
for domain = domains;
    filestring = dir([partype num2str(domain) '.mat']);
    sizevar = size(filestring);
    if sizevar(1) == 0
        loadps = 1;
    else
        loadps = 0;
    end
    if loadps
        filestring = dir([partype '*.*.' num2str(domain)]);
        filestring = filestring.name;
        [ ~, px, ~, pz, pux, puy, puz, pq ] =...
            load_domain_particles([directory '\'  filestring]);
        save([partype num2str(domain) '.mat'],...
            'px','pz','pux','puy','puz','pq');
    else
        load([partype num2str(domain) '.mat']);
    end
    mx=mean(px);

    pbx = interp2(xv,zv,bx',px,pz);
    pby = interp2(xv,zv,by',px,pz);
    pbz = interp2(xv,zv,bz',px,pz);



    Bmod=sqrt(pbx.^2+pby.^2+pbz.^2);
    bxu=pbx./Bmod;
    byu=pby./Bmod;
    bzu=pbz./Bmod;


    bp1x=bzu;
    bp1z=-bxu; 
    Bmodp=sqrt(bp1x.^2+bp1z.^2);
    bp1x=bp1x./Bmodp;
    bp1z=bp1z./Bmodp;
    bp1y=bp1z*0;


    %make  unit vector perp to B and bp1 
    bp2x=byu.*bp1z-bzu.*bp1y;
    bp2y=bzu.*bp1x-bxu.*bp1z;
    bp2z=bxu.*bp1y-byu.*bp1x;



    Vpar = bxu.*pux+byu.*puy+bzu.*puz;
    %Vpoff=[-0.2743   -0.0533];
    Vperp1 = bp1x.*pux+bp1y.*puy+bp1z.*puz;
    Vperp2 = bp2x.*pux+bp2y.*puy+bp2z.*puz;

    Vperp=sqrt(Vperp1.^2+Vperp2.^2);


    cntx=cntx+1;
    Xps(cntx) = mean(px);

    cntz=0;
    for z0 = zs
        cntz=cntz+1;

        dz = 1; %potentially important, work into varargin

        q = pq(abs(pz-z0)<dz); %index
        ux = Vpar(abs(pz-z0)<dz);
        uy = Vperp(abs(pz-z0)<dz);
        
        up1 = Vperp1(abs(pz-z0)<dz);
        up2 = Vperp2(abs(pz-z0)<dz);


        Nvx=50;
        Nvy=51;

        vxs = linspace(-Vlim,Vlim,Nvx); %vpar
        vys = linspace(0,Vlim,Nvy); %%%%%%%%%%%%% %vperp
        
        vp1s = linspace(-Vlim,Vlim,Nvx);

        [VX,VY]= ndgrid(vxs,vys);
        [VP1,VP2] = ndgrid(vp1s,vp1s);


        dvx = vxs(2)-vxs(1);
        dvy = vys(2)-vys(1);

        xsub = round((ux-min(vxs))/dvx)+1; xsub = min(xsub,Nvx); xsub = max(xsub,ones(size(xsub)));
        ysub = round((uy-min(vys))/dvy)+1; ysub = min(ysub,Nvx); ysub = max(ysub,ones(size(xsub)));
        
        p1sub = round((up1-min(vp1s))/dvx)+1; p1sub = min(p1sub,Nvx); p1sub = max(p1sub,ones(size(p1sub)));
        p2sub = round((up2-min(vp1s))/dvx)+1; p2sub = min(p2sub,Nvx); p2sub = max(p2sub,ones(size(p2sub)));

        Fxy = accumarray([xsub',ysub'],-q,[Nvx Nvy])./abs(VY);
        Fpp = accumarray([p1sub',p2sub'],-q,[Nvx Nvx]);

        %colbnds = [0 2]*1e-3;
        %colbnds=[-4.5 -1]-1;
        % figure(300),subplot(length(zs),length(doms),length(doms)*length(zs)-(cntz-1)*length(doms)-cntx+1)
        % pcolor(vxs,vys,log10(Fxyt')),shading flat
        % caxis(colbnds)
        % axis off
        % colormap(jet)
        % xlabel('vx')
        % ylabel('vy')

        % subplot(2,3,2)
        % pcolor(vxs,vzs,Fxzt'),shading flat
        % caxis(colbnds)
        % xlabel('vx')
        % ylabel('vz')
        % title('Sphere Electrons')
        % 
        % subplot(2,3,3)
        % pcolor(vys,vzs,Fyzt'),shading flat
        % caxis(colbnds)
        % xlabel('vy')
        % ylabel('vz')

        %
        %colbnds = [0 2]*1e-2;
        %colbnds=[-4.5 -1];
        %figure(301),subplot(length(zs)/2,length(doms)*2,length(doms)*length(zs)-(cntz-1)*length(doms)-cntx+1),cla

        %axes(Hvec(cntz))
        figure(floor((cntz-1)/9)+1)
        subplot(3,3,mod(cntz-1,9)+1)
        hold off
        pcolor(vxs,vys,log10(Fxy)'),shading interp
        xlabel('v_{||}');ylabel('v_{\perp}');
        title(num2str(z0));
        hold on
        caxis(colbnds)
        axis off
        axis([-Vlim Vlim 0 Vlim])
%         colormap(hot)
        hold on

        contour(vxs,vys,sqrt(VX.^2+VY.^2)',[-0.5 0.3 .4 .5],'m')

        plot([0 0],[0 1],'m')


        figure(floor((cntz-1)/9)+101)
        subplot(3,3,mod(cntz-1,9)+1)
        hold off
        pcolor(vp1s,vp1s,log10(Fpp)'),shading interp
        xlabel('v_{\perp 1}');ylabel('v_{\perp 2}');
        title(num2str(z0));
        hold on
        caxis(colbndsperp)
        axis off
        axis([-Vlim Vlim -Vlim Vlim])
%         colormap(hot)
        hold on

        contour(vp1s,vp1s,sqrt(VP1.^2+VP2.^2)',[-0.5 0.3 .4 .5],'m')

        plot([-1 1],[0 0],'m')
        plot([0 0],[-1 1],'m')
        % pause(0.1)
    end


end

% title(num2str(mx/20))

%axes(H2)
%
figure()
mpx =mean(px);
plot(mpx/sqrt(mime),zs/sqrt(mime),'xk','markersize',4,'linewidth',3)

%%
% figure(30),clf
% pcolor(xv,zv,ez'),shading flat, hold on,colorbar
% caxis([-0.1 0.1])
% % caxis([-0.2 0.1])
% % caxis([-1 3])
% axis([-40 100 -40 40])
% contour(xv,zv(iz:end),Ay(:,iz:end)',[-zs-0.1 -zs+0.1],'m')
% % contour(xv,zv,emix',0.99*[-1 1],'k')
% title('E_z')
% xlabel('x/de')
% ylabel('z/de')
% colormap(jet)
% 
% for jj=1:length(Xps)
%     plot([-0.5 -0.5]+Xps(jj)-200,[-40 40],'k')
%     plot([0.5 0.5]+Xps(jj)-200,[-40 40],'k')
%    
% 
% end
% 
