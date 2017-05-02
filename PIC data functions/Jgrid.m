function Jvals = Jgrid(xvals, zvals, vvec,timestep)

    load(['Ex_' num2str(timestep) '.mat']);
    load(['Ey_' num2str(timestep) '.mat']);
    load(['Ez_' num2str(timestep) '.mat']);
    load(['Bx_' num2str(timestep) '.mat']);
    load(['By_' num2str(timestep) '.mat']);
    load(['Bz_' num2str(timestep) '.mat']);
    load(['magB_' num2str(timestep) '.mat']);
    load(['Psi_' num2str(timestep) '.mat']);
    
    earrs = {ex, ey, ez};
    barrs = {bx,by,bz};
    
    nx = length(xvals);
    nz = length(zvals);
    
    Jvals = zeros(nx,nz);
    
    for q = 1:(nx*nz)
        [i,j] = ind2sub([nx,nz],q);
        Jvals(i,j) = invariantJCmat2([xvals(i) zvals(j) vvec],earrs,barrs,magb,Psi);
    end
    
%     for i = 1:nx
%         for j = 1:nz
%             Jvals(i,j) = invariantJCmat2([xvals(i) zvals(j) vvec],earrs,barrs,magb,Psi);
%         end
%     end
            
end    