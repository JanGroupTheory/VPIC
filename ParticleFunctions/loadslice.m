function qv = loadslice(q,slice,nx,nz,datadir)


switch q
    case 'Pexx'; q='pe-xx';
    case 'Pexy'; q='pe-xy';
    case 'Pexz'; q='pe-xz';
    case 'Peyy'; q='pe-yy';
    case 'Peyz'; q='pe-yz';
    case 'Pezz'; q='pe-zz';
    case 'Pixx'; q='pi-xx';
    case 'Pixy'; q='pi-xy';
    case 'Pixz'; q='pi-xz';
    case 'Piyy'; q='pi-yy';
    case 'Piyz'; q='pi-yz';
    case 'Pizz'; q='pi-zz';
    case 'emix'; q= 'e-mix1';
end
   
   fd = fopen([ datadir '\' q '.gda'],['rb', ]);
   for nslice = 1:slice
       qv = fread(fd,[nx nz],'single')';
   end
   
   %decimate
   qv = (qv(1:2:end,:) + qv(2:2:end,:))/2;
   qv = (qv(:,1:2:end) + qv(:,2:2:end))/2;

   
   fclose(fd);