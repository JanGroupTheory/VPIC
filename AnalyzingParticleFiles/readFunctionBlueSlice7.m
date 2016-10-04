%Reads in PIC data from '\\VTF9\share\Daughton\Jan3\'.
%Options for quanitities are qs = {'Ay','Bx','By','Bz','Ex','Ey','Ez','Jx','Jy','Jz','ne',...
%    'Pe-xx','Pe-xy','Pe-xz','Pe-yy','Pe-yz','Pe-zz','rho','Uex','Uey','Uez','Uix','Uiy','Uiz'}
%Choose slice 0 to 215


function [temp1] = readFunctionBlueSlice7(tslice,name,savedir)

readdir=savedir;

str = sprintf('%1$s_%2$s.mat', char(name),num2str(tslice));
tmp = load(sprintf('%1$s%2$s', readdir,str));

if strcmp(char(name),'Pepar')
    temp1 = tmp.Pepar;
elseif strcmp(char(name), 'Peperp1')
    temp1 = tmp.Peperp1;
elseif strcmp(char(name), 'Peperp2')
    temp1 = tmp.Peperp2;
elseif strcmp(char(name),'Pipar')
    temp1 = tmp.Pipar;
elseif strcmp(char(name), 'Piperp1')
    temp1 = tmp.Piperp1;
elseif strcmp(char(name), 'Piperp2')
    temp1 = tmp.Piperp2;
elseif strcmp(char(name), 'Ay')
    temp1 = tmp.ay;
elseif strcmp(char(name), 'Phi')
    temp1 = tmp.phi;
elseif strcmp(char(name), 'Psi')
    temp1 = tmp.Psi;
elseif strcmp(char(name), 'Phipa')
    temp1 = tmp.Phipa;
elseif strcmp(char(name), 'Ppar')
    temp1 = tmp.Ppar;
elseif strcmp(char(name), 'Pperp1')
    temp1 = tmp.Pperp1;
elseif strcmp(char(name), 'Pperp2')
    temp1 = tmp.Pperp2;
elseif strcmp(char(name), 'Pparp1')
    temp1 = tmp.Pparp1;
elseif strcmp(char(name), 'Pparp2')
    temp1 = tmp.Pparp2;
elseif strcmp(char(name), 'Pp1p2')
    temp1 = tmp.Pp1p2;
elseif strcmp(char(name), 'Buppoints')
    temp1 = tmp.buppoints;
elseif strcmp(char(name), 'Bdownpoints')
    temp1 = tmp.bdownpoints;
elseif strcmp(char(name), 'Bx')
    temp1 = tmp.bx;
elseif strcmp(char(name), 'By')
    temp1 = tmp.by;
elseif strcmp(char(name), 'Bz')
    temp1 = tmp.bz;
elseif strcmp(char(name), 'Jx')
    temp1 = tmp.jx;
elseif strcmp(char(name), 'Jy')
    temp1 = tmp.jy;
elseif strcmp(char(name), 'Jz')
    temp1 = tmp.jz;
elseif strcmp(char(name), 'Ex')
    temp1 = tmp.ex;
elseif strcmp(char(name), 'Ey')
    temp1 = tmp.ey;
elseif strcmp(char(name), 'Ez')
    temp1 = tmp.ez;
elseif strcmp(char(name), 'Uex')
    temp1 = tmp.uex;
elseif strcmp(char(name), 'Uey')
    temp1 = tmp.uey;
elseif strcmp(char(name), 'Uez')
    temp1 = tmp.uez;
elseif strcmp(char(name), 'Uix')
    temp1 = tmp.uix;
elseif strcmp(char(name), 'Uiy')
    temp1 = tmp.uiy;
elseif strcmp(char(name), 'Uiz')
    temp1 = tmp.uiz;
elseif strcmp(char(name), 'Ne')
    temp1 = tmp.ne;
elseif strcmp(char(name), 'Ni')
    temp1 = tmp.ni;
elseif strcmp(char(name), 'PtenTotal')
    temp1 = tmp.PtenTotal;
elseif strcmp(char(name), 'Pxx')
    temp1 = tmp.Pxx;
elseif strcmp(char(name), 'Pyy')
    temp1 = tmp.Pyy;
elseif strcmp(char(name), 'Pzz')
    temp1 = tmp.Pzz;
elseif strcmp(char(name), 'Pxy')
    temp1 = tmp.Pxy;
elseif strcmp(char(name), 'Pxz')
    temp1 = tmp.Pxz;
elseif strcmp(char(name), 'Pyz')
    temp1 = tmp.Pyz;
else
    temp1 = tmp.temptemp;
end

end






