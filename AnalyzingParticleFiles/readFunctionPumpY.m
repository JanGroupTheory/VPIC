%Options for quanitities are qs = {'Ay','Bx','By','Bz','Ex','Ey','Ez','Jx','Jy','Jz','ne',...
%    'Pe-xx','Pe-xy','Pe-xz','Pe-yy','Pe-yz','Pe-zz','rho','Uex','Uey','Uez','Uix','Uiy','Uiz'}
%Choose slice 0 to 95

function [temp] = readFunctionPumpY(slice, name, readdir, nlx, nlz)

%data followed by time and some other tag
startbyte = (slice-1)*4*(nlx*nlz+0);
% startbyte = (slice-1)*4*(nlx*nlz+2);
% startbyte = (slice-1)*4*(nlx*nlz + 0);

filenam=[char(readdir)  char(name) '.gda'];

%filenam

fid = fopen(filenam,'r');
fseek(fid,startbyte,'bof');
temp = fread(fid,[nlx,nlz],'float32');
fclose(fid);

%decimate
temp = (temp(1:2:end,:) + temp(2:2:end,:) )/2;
temp = (temp(:,1:2:end) + temp(:,2:2:end) )/2;



