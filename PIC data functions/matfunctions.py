import numpy as np
import scipy.io as sio
from io import StringIO
import os.path
import os
from os import listdir, mkdir
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import scipy as sp
import scipy.ndimage
import scipy.interpolate
import multiprocessing
from itertools import repeat
import sys
from matplotlib.animation import FuncAnimation
#from images2gif import writeGIF
from PIL import Image
#import moviepy.editor as mpy
def saveslice(var,slicenumP,slicenumM,nx,nz,gdadir,savedir):
    casedict = {
        'pe-xx':'Pexx',
        'pe-xy':'Pexy',
        'pe-xz':'Pexz',
        'pe-yy':'Peyy',
        'pe-yz':'Peyz',
        'pe-zz':'Pezz',
        'pi-xx':'Pixx',
        'pi-xy':'Pixy',
        'pi-xz':'Pixz',
        'pi-yy':'Piyy',
        'pi-yz':'Piyz',
        'pi-zz':'Pizz'
    }
    var = casedict.get(var,var)
    outfile = savedir + '/' + var +'_' + str(slicenumM)  +'.mat'
    #outfile  = savedir + '/' + var + str(slicenumM) + '.mat'

    if not (os.path.isfile(outfile)):
        out = loadslice(var,slicenumP,nx,nz,gdadir)
        out = out[::2,:]+out[1:out.shape[0]:2,:] #decimation steps
        out = out[:,::2]+out[:,1:out.shape[1]:2]
        sio.savemat(outfile,{var:out/4})
    return 0

def loadslice(q,slicenum,nx,nz,datadir):
    q = ''.join(q.split('_ave'))
    casedict = {
        'Pexx':'pe-xx',
        'Pexy':'pe-xy',
        'Pexz':'pe-xz',
        'Peyy':'pe-yy',
        'Peyz':'pe-yz',
        'Pezz':'pe-zz',
        'Pixx':'pi-xx',
        'Pixy':'pi-xy',
        'Pixz':'pi-xz',
        'Piyy':'pi-yy',
        'Piyz':'pi-yz',
        'Pizz':'pi-zz',
        'neUUexx':'neuue-xx',
        'neUUexy':'neuue-xy',
        'neUUexz':'neuue-xz',
        'neUUeyy':'neuue-yy',
        'neUUeyz':'neuue-yz',
        'neUUezz':'neuue-zz'
    }
    q = casedict.get(q,q)
    
    with open(datadir+'/'+q+'_'+str(slicenum)+'.gda','rb') as fd:
        return np.reshape(np.fromfile(fd,dtype='f4',count = nx*nz),(nz,nx))

def process1slice(matdir,slicenum,dt,dx_de):
    mfcont = sio.loadmat(matdir+'/Slice_'+str(slicenum)+'.mat')
    if 'ne' in mfcont:
        ne = mfcont['ne']
        ni = mfcont['ni']
        bx = mfcont['bx']
        by = mfcont['by']
        bz = mfcont['bz']
        ey = mfcont['ey']
        absB = mfcont['absB']
        Pexx = mfcont['Pexx']
        Pexy = mfcont['Pexy']
        Pexz = mfcont['Pexz']
        Peyy = mfcont['Peyy']
        Peyz = mfcont['Peyz']
        Pezz = mfcont['Pezz']
        Pixx = mfcont['Pixx']
        Pixy = mfcont['Pixy']
        Pixz = mfcont['Pixz']
        Piyy = mfcont['Piyy']
        Piyz = mfcont['Piyz']
        Pizz = mfcont['Pizz']
        bhatx = np.divide(bx,absB)
        bhaty = np.divide(by,absB)
        bhatz = np.divide(bz,absB)
        Ppar = Pexx*bhatx**2+Peyy*bhaty**2+Pezz*bhatz**2 + 2*(bhatx*(Pexy*bhaty+Pexz*bhatz)+bhatz*bhaty*Peyz)
        Pipar = Pixx*bhatx**2+Piyy*bhaty**2+Pizz*bhatz**2 + 2*(bhatx*(Pixy*bhaty+Pixz*bhatz)+bhatz*bhaty*Piyz)
        ehat1x = -bhaty/np.sqrt(1-bhatz**2)
        ehat1y = bhatx/np.sqrt(1-bhatz**2)
        ehat1z = 0*bhaty
        ehat2x = -bhatz*ehat1y 
        ehat2y = bhatz*ehat1x 
        ehat2z = bhatx*ehat1y-bhaty*ehat1x
        a = Pexx*ehat1x**2+Peyy*ehat1y**2+Pezz*ehat1z**2 + 2*(ehat1x*(Pexy*ehat1y+Pexz*ehat1z)+ehat1z*ehat1y*Peyz);
        d = Pexx*ehat2x**2+Peyy*ehat2y**2+Pezz*ehat2z**2 + 2*(ehat2x*(Pexy*ehat2y+Pexz*ehat2z)+ehat2z*ehat2y*Peyz);
        b = Pexx*ehat1x*ehat2x+Peyy*ehat1y*ehat2y+Pezz*ehat1z*ehat2z + Pexy*(ehat1x*ehat2y+ehat2x*ehat1y) \
            + Pexz*(ehat1x*ehat2z+ehat2x*ehat1z)+ Peyz*(ehat1z*ehat2y+ehat2z*ehat1y)
        Pperp1 = (a+d)/2+np.sqrt(((a-d)/2)**2+b**2) 
        Pperp2 = (a+d)/2-np.sqrt(((a-d)/2)**2+b**2) 
        Pp1p2 = 0*Ppar
        p1 = Pexx*ehat1x*bhatx+Peyy*ehat1y*bhaty+Pezz*ehat1z*bhatz + Pexy*(ehat1x*bhaty+bhatx*ehat1y) \
             + Pexz*(ehat1x*bhatz+bhatx*ehat1z)+ Peyz*(ehat1z*bhaty+bhatz*ehat1y);
        p2 = Pexx*ehat2x*bhatx+Peyy*ehat2y*bhaty+Pezz*ehat2z*bhatz + Pexy*(ehat2x*bhaty+bhatx*ehat2y) \
             + Pexz*(ehat2x*bhatz+bhatx*ehat2z)+ Peyz*(ehat2z*bhaty+bhatz*ehat2y);
        theta = .5*np.arcsin(2*b/(a+d))
        Pparp1 = p1*np.cos(theta)-p2*np.sin(theta)
        Pparp2 = p2*np.cos(theta)+p1*np.sin(theta)
        
        a = Pixx*ehat1x**2+Piyy*ehat1y**2+Pizz*ehat1z**2 + 2*(ehat1x*(Pixy*ehat1y+Pixz*ehat1z)+ehat1z*ehat1y*Piyz)
        d = Pixx*ehat2x**2+Piyy*ehat2y**2+Pizz*ehat2z**2 + 2*(ehat2x*(Pixy*ehat2y+Pixz*ehat2z)+ehat2z*ehat2y*Piyz)
        b = Pixx*ehat1x*ehat2x+Piyy*ehat1y*ehat2y+Pizz*ehat1z*ehat2z + Pixy*(ehat1x*ehat2y+ehat2x*ehat1y) \
            + Pixz*(ehat1x*ehat2z+ehat2x*ehat1z)+ Piyz*(ehat1z*ehat2y+ehat2z*ehat1y)
        Piperp1 = (a+d)/2+np.sqrt(((a-d)/2)**2+b**2) 
        Piperp2 = (a+d)/2-np.sqrt(((a-d)/2)**2+b**2) 
        Pip1p2 = 0*Ppar
        p1 = Pixx*ehat1x*bhatx+Piyy*ehat1y*bhaty+Pizz*ehat1z*bhatz + Pixy*(ehat1x*bhaty+bhatx*ehat1y) \
             + Pixz*(ehat1x*bhatz+bhatx*ehat1z)+ Piyz*(ehat1z*bhaty+bhatz*ehat1y)
        p2 = Pixx*ehat2x*bhatx+Piyy*ehat2y*bhaty+Pizz*ehat2z*bhatz + Pixy*(ehat2x*bhaty+bhatx*ehat2y) \
             + Pixz*(ehat2x*bhatz+bhatx*ehat2z)+ Piyz*(ehat2z*bhaty+bhatz*ehat2y)
        theta = .5*np.arcsin(2*b/(a+d))
        Piparp1 = p1*np.cos(theta)-p2*np.sin(theta)
        Piparp2 = p2*np.cos(theta)+p1*np.sin(theta)
        
        Tpar = Ppar/ne
        Tperp = (Pperp1+Pperp2)/2/ne
        Tipar = Pipar/ni 
        Tiperp = (Piperp1+Piperp2)/2/ni
        
        #eycorner = np.mean(np.mean(ey[1:10,1:10])))

        sz = np.shape(bx)
        Ibx = np.cumsum(bz[0,:])-(bz[0,:]+bz[0,0])/2
        Ibz = np.cumsum(bx, axis=0) - (np.matmul(np.ones((sz[0],1)),np.reshape(bx[0,:],(1,-1))) + bx)/2
        Psi = 2*(np.matmul(np.ones((sz[0],1)), np.reshape(Ibx,(1,-1))) - Ibz)*dx_de

        #Psi = (-np.cumsum(bx,axis=0)+np.cumsum(bz,axis=1))*dx_de;
        mfcont.update({'Ppar':Ppar,'Pperp1':Pperp1,'Pperp2':Pperp2,'Pparp1':Pparp1,
                       'Pparp2':Pparp2,'Pp1p2':Pp1p2,'Pipar':Pipar,'Piperp1':Piperp1,
                       'Piperp2':Piperp2,'Piparp1':Piparp1,'Piparp2':Piparp2,'Pip1p2':Pip1p2,
                       'Psi':Psi,'Tpar':Tpar,'Tperp':Tperp,'Tipar':Tipar,'Tiperp':Tiperp})
        
        sio.savemat(matdir+ '/Slice_'+ str(slicenum)+ '.mat',mfcont)
    else:
        print('Slice ' + str(slicenum) + ' is missing non-averaged quantities!\n')
def processSlice(matdir, Slicenums):
    with open(matdir + '/../../info','r') as fp:
        info = fp.read()

    mime = int(float(info.split('mi/me =')[1].split()[0]))
    dx_de = float(info.split('dx/de =')[1].split()[0])
    wpewce = float(info.split('wpe/wce =')[1].split()[0])
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    #v_A = float(info.split('v_A')[1].split('=')[1].split()[0])
    dt0 = wpewce*mime
    dt = int(1/dtwci)
    Edrive = float(info.split('edrive =')[1].split()[0])
    taudrive = float(info.split('tdrive =')[1].split()[0])

    #taudrive = dt0*tauwci
    #Edrive = v_A/wpewce*edrivefactor
    #dt = mime*wpewce

    #for slicenum in Slicenums:
     #   process1slice(matdir,slicenum,dt,dx_de)
    pool = multiprocessing.Pool()
    pool.starmap(process1slice,zip(repeat(matdir),Slicenums,repeat(dt),repeat(dx_de)))
    fixPsi(matdir,Slicenums,dt0, Edrive, taudrive)
    return 0

def fixPsi(matdir, Slicenums, dt, Edrive, taudrive):
    Psi0 = 0;
    for slicenum in Slicenums:
         mfcont = sio.loadmat(matdir+'/Slice_'+str(slicenum)+'.mat')
         #EY = mfcont['eycorner'];
         #Psi0 -= EY*dt
         t0 = dt*slicenum
         Psi0 -= Edrive*(dt+taudrive*(1-np.exp(dt/taudrive))*np.exp(-t0/taudrive))
         mfcont.update({'Psi':mfcont['Psi']+Psi0,'Psi0':Psi0})
         sio.savemat(matdir+'/Slice_'+str(slicenum)+'.mat',mfcont)

def fixPsi_ave(matdir, slicenum):
    mfcont = sio.loadmat(matdir+'/Slice_'+str(slicenum)+'.mat')
    if 'Psi_ave' in mfcont:
        mfcont.update({'Psi_ave':mfcont['Psi_ave']+mfcont['Psi0']})
        sio.savemat(matdir+'/Slice_'+str(slicenum)+'.mat',mfcont)
    else:
        print('Slice ' + str(slicenum) + 'is missing Psi_ave!')

def bundleslices(basedir,slicenums):
    indir = basedir + '/Slices/'
    gdadir = basedir + '/data/'

    if not os.path.isdir(indir+'FullSlice/'):
        mkdir(indir + 'FullSlice/')
    
    varlist = []
    varz = [f for f in listdir(indir) if ((os.path.isfile(indir+f)) and (f!='info'))]
    varz = np.sort(varz)
    length = len(varz)
    ts = []
    i = 0
    flag = True
    while flag:
        temp = varz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(ts))==0:
            ts.append(tnew)
        else:
            flag = False
    ts = np.sort(ts).astype(int)
    j = -1
    for i in range(0,int(length/len(ts))):
        temp = varz[i*len(ts)].split('_')
        varlist.append(temp[0])

    #for slicenum in slicenums:
     #   bundler(indir,varlist,slicenum)
    pool = multiprocessing.Pool()
    pool.starmap(bundler,zip(repeat(indir),repeat(varlist),slicenums))
    return 0

def bundler(indir,varlist,slicenum):
    savedict = {}
    for i in range(0,len(varlist)):
        savedict.update(sio.loadmat(indir+ '/'+ varlist[i]+ '_'+ str(slicenum)+ '.mat'))
        #savedict.update(sio.loadmat(indir+ '/' + varlist[i]+ str(slicenum)+ '.mat'))
    outfile = indir + 'FullSlice/Slice_'+ str(slicenum) + '.mat'
    if not (os.path.isfile(outfile)):
        sio.savemat(outfile, savedict)

def populateslices(varlist,ts,nx,nz,gdadir,savedir):
    pool = multiprocessing.Pool()
    for var in varlist:
        pool.starmap(saveslice, zip(repeat(var),ts,range(0,len(ts)),repeat(nx),repeat(nz),repeat(gdadir),repeat(savedir)))
        #for slicenumM in range(0,len(ts)-1):
         #   saveslice(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir)
        #Parallel(n_jobs=-1)(delayed(saveslice)(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir) for slicenumM in range(0,len(ts)-1)) #now parallelized over time slice
    return 0

#def populateslicesvarpar(varlist,ts,nx,nz,gdadir,savedir):
 #   for slicenumM in range(0,len(ts)-1):
  #      Parallel(n_jobs=-1)(delayed(saveslice)(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir) for var in varlist) #now parallelized over variable for smaller number of time slices
   # return 0

def calcalpha12(indir,t0,t1,nPsi,vA,wpewce,mime):
    m = 1
    e = -1
    dtwpe = wpewce*mime
    alpha1 = 0.0
    alpha2 = 0.0
    slicevals = sio.loadmat(indir+'/Slice_'+str(t0)+'.mat')
    Jmat0 = slicevals['Jmat']
    mumat0 = slicevals['mumat']
    slicevals = sio.loadmat(indir+'/Slice_'+str(t1)+'.mat')
    Jmat1 = slicevals['Jmat']
    mumat1 = slicevals['mumat']
    nU,nmu,_ = np.shape(Jmat0)
    mumax = np.zeros((nU,1))
    mumin = np.zeros((nU,1))
    J0 = np.zeros(np.shape(Jmat0[:,:,0]))
    J1 = np.zeros(np.shape(Jmat1[:,:,0]))
    for j in range(0,np.size(nPsi)):
        for i in range(0,nU):
            k = nPsi[j]
            mumax = min((mumat0[i,nmu-1,k],mumat1[i,nmu-1,k]))
            mumin = max((mumat0[i,0,k],mumat1[i,0,k]))
            mu = np.linspace(mumin,mumax,nmu)
            J0[i,:] = np.interp(mu,mumat0[i,:,k],Jmat0[i,:,k])
            J1[i,:] = np.interp(mu,mumat1[i,:,k],Jmat1[i,:,k])
        dJdt = (J1-J0)/((t1-t0)*dtwpe)
        alpha1 += -np.sqrt(m/2)/vA*np.mean(np.mean(dJdt,1)/np.sqrt(slicevals['U'].T))
        alpha2 += np.sqrt(m/2/vA**2*np.mean(np.mean(dJdt**2,1)/slicevals['U'].T))
    alpha1 /= np.size(nPsi)
    alpha2 /= np.size(nPsi)
    return alpha1, alpha2

def reinterpJ(Jmat0,mumat0,Jmat1,mumat1):
    return J0, J1

def calcJ(indir,slicenums,Umax,nU,nmu):
    m = 1 
    e = -1
    Us = np.linspace(0,Umax,num=nU+1)
    Us = Us[1:np.size(Us)]
    for slicenum in slicenums:
        slicevals=sio.loadmat(indir+ '/Slice_'+str(slicenum)+ '.mat')
        nPsi = np.size(slicevals['Psivals'])
        Jmat = np.zeros((nU,nmu,nPsi))
        mumat = np.zeros((nU,nmu,nPsi))
        for k in range(0,nPsi-1):
            l = slicevals['l'+str(k)].flatten()
            dl = (np.append(0,np.diff(l))+np.append(np.diff(l),0))/2
            phipar = slicevals['phipar'+str(k)].flatten()
            B = slicevals['B'+str(k)].flatten()
            #Bm = np.interp(0,l,B)
            Bm = np.min(B)
            Binf = np.max(B)
            phiparm = -phipar[0]
            for i in range(0,nU):
                muc = (Us[i]+e*phiparm)/Binf
                if muc<0:
                    muc = 0
                mum = Us[i]/Bm
                mu = np.linspace(muc,mum,nmu)
                mumat[i,:,k] = mu
                for j in range(0,nmu):
                    Enerpar = np.nan_to_num(np.sqrt(2/m*(Us[i]-mu[j]*B-e*phipar)))
                    Jmat[i,j,k] = 2*np.sum(Enerpar*np.abs(dl))
        U = Us
        dic = sio.loadmat(indir + '/Slice_'+ str(slicenum)+ '.mat')
        dic.update({'Jmat':Jmat,'U':U,'mumat':mumat})
        sio.savemat(indir + '/Slice_'+ str(slicenum)+ '.mat', dic)

def makemats(basedir):
    savedir = basedir + '/Slices/'
    gdadir = basedir + '/data/'

    if not os.path.isdir(savedir):
        mkdir(savedir)
    
    varlist = []
    varz = [f for f in listdir(gdadir) if ((os.path.isfile(gdadir+f)) and (f!='info'))]
    varz = np.sort(varz)
    length = len(varz)
    ts = []
    i = 0
    flag = True
    while flag:
        temp = varz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(ts))==0:
           ts.append(tnew)
        else:
            flag = False
    ts = np.sort(ts).astype(int)
    j = -1
    for i in range(0,int(length/len(ts))):
        temp = varz[i*len(ts)].split('_')
        varlist.append(temp[0])
    with open(basedir + '/info','r') as fp:
        info = fp.read()

    nsteps = int(info.split('num_step =')[1].split()[0])
    nx = int(float(info.split('nx =')[1].split()[0]))
    nz = int(float(info.split('nz =')[1].split()[0]))
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    nslices = int(nsteps*dtwci);
   # nx = 2560;
   # nz = 2560;
    fp.close()
    populateslices(varlist,ts,nx,nz,gdadir,savedir)
    bundleslices(basedir,range(0,len(ts)))
    processSlice(savedir+'/FullSlice/',range(0,len(ts)))

def contourvalues(indir, slicenums, xv, zv, Psivals, mime):
    ns = 10
    for slicenum in slicenums:
        fdict = sio.loadmat(indir + '/Slice_' + str(slicenum) + '.mat')
        absB = fdict['absB']
        bhatx = fdict['bx']/absB
        bhaty = fdict['by']/absB
        bhatz = fdict['bz']/absB

        ex = sp.ndimage.filters.gaussian_filter(fdict['ex'],ns,mode='constant')
        ey = sp.ndimage.filters.gaussian_filter(fdict['ey'],ns,mode='constant')
        ez = sp.ndimage.filters.gaussian_filter(fdict['ez'],ns,mode='constant')

        Epar = ex*bhatx+ey*bhaty+ez*bhatz

        splineB = sp.interpolate.RectBivariateSpline(xv,zv,absB.T)
        splineE = sp.interpolate.RectBivariateSpline(xv,zv,Epar.T)
        splineBg = sp.interpolate.RectBivariateSpline(xv,zv,bhaty.T)
        splinebx = sp.interpolate.RectBivariateSpline(xv,zv,bhatx.T)
        splinebz = sp.interpolate.RectBivariateSpline(xv,zv,bhatz.T)
        newdict = {}
        for j in range(0,len(Psivals)): 
            temp = plt.contour(xv,zv,fdict['Psi'],Psivals).collections[0].get_paths()[0].vertices
            x = temp[:,0]
            z = temp[:,1]
            B = splineB.ev(x,z)
            B[B!=B] = 1
        
            epar = splineE.ev(x,z)
            epar[epar!=epar]  = 0
            dx = (np.append(0,np.diff(x))+np.append(np.diff(x),0))/2
            dz = (np.append(0,np.diff(z))+np.append(np.diff(z),0))/2

            bg = splineBg.ev(x,z)
            bg[bg!=bg] = 0
            bx = splinebx.ev(x[0],z[0])
            bz = splinebz.ev(x[0],z[0])
            parsign = np.sign(bx*dx[0]+bz*dz[0])
            dl = parsign*np.sqrt((dx**2+dz**2)*(1/(1-bg**2)))*np.sqrt(mime)
            l = np.cumsum(dl)
            phipar = -np.cumsum(dl*epar)
            midind = np.argmin(np.abs(z))
            l = l - l[midind]
            phipar = phipar - phipar[midind]
            newdict.update({'x'+str(j):x,'z'+str(j):z,'l'+str(j):l,'B'+str(j):B,'phipar'+str(j):phipar})
            dic = sio.loadmat(indir + '/Slice_'+ str(slicenum)+ '.mat') 
            dic.update({'xv':xv,'zv':zv,'Psivals':Psivals})
            dic.update(newdict)
            sio.savemat(indir + '/Slice_'+ str(slicenum)+ '.mat', dic)

def saveavgslice(var,slicenumP,slicenumM,nx,nz,gdadir,savedir):
    casedict = {
        'pe-xx':'Pexx',
        'pe-xy':'Pexy',
        'pe-xz':'Pexz',
        'pe-yy':'Peyy',
        'pe-yz':'Peyz',
        'pe-zz':'Pezz',
        'pi-xx':'Pixx',
        'pi-xy':'Pixy',
        'pi-xz':'Pixz',
        'pi-yy':'Piyy',
        'pi-yz':'Piyz',
        'pi-zz':'Pizz',
        'neuue-xx':'neUUexx',
        'neuue-xy':'neUUexy',
        'neuue-xz':'neUUexz',
        'neuue-yy':'neUUeyy',
        'neuue-yz':'neUUeyz',
        'neuue-zz':'neUUezz'
    }
    var = casedict.get(var,var)+'_ave'
    outfile = savedir + '/' + var +'_' + str(slicenumM)  +'.mat'
    #outfile  = savedir + '/' + var + str(slicenumM) + '.mat'

    if not (os.path.isfile(outfile)):
        out = loadslice(var,slicenumP,nx,nz,gdadir) #update dictionary in loadslice
        out = out[::2,:]+out[1:out.shape[0]:2,:] #decimation steps
        out = out[:,::2]+out[:,1:out.shape[1]:2]
        sio.savemat(outfile,{var:out/4})
    return 0

def process1sliceavg(matdir,slicenum,dt,dx_de):
    mfcont = sio.loadmat(matdir+'/Slice_'+str(slicenum)+'.mat')
    if 'ne_ave' in mfcont:
        ne = mfcont['ne_ave']
        bx = mfcont['bx_ave']
        by = mfcont['by_ave']
        bz = mfcont['bz_ave']
        ey = mfcont['ey_ave']
        absB = np.sqrt(bx**2+by**2+bz**2)
        Pexx = mfcont['Pexx_ave']
        Pexy = mfcont['Pexy_ave']
        Pexz = mfcont['Pexz_ave']
        Peyy = mfcont['Peyy_ave']
        Peyz = mfcont['Peyz_ave']
        Pezz = mfcont['Pezz_ave']
        neUUexx = mfcont['neUUexx_ave']
        neUUexy = mfcont['neUUexy_ave']
        neUUexz = mfcont['neUUexz_ave']
        neUUeyy = mfcont['neUUeyy_ave']
        neUUeyz = mfcont['neUUeyz_ave']
        neUUezz = mfcont['neUUezz_ave']
        bhatx = np.divide(bx,absB)
        bhaty = np.divide(by,absB)
        bhatz = np.divide(bz,absB)
        Ppar = Pexx*bhatx**2+Peyy*bhaty**2+Pezz*bhatz**2 + 2*(bhatx*(Pexy*bhaty+Pexz*bhatz)+bhatz*bhaty*Peyz)
        neUUepar = neUUexx*bhatx**2+neUUeyy*bhaty**2+neUUezz*bhatz**2 + 2*(bhatx*(neUUexy*bhaty+neUUexz*bhatz)+bhatz*bhaty*neUUeyz)
        ehat1x = -bhaty/np.sqrt(1-bhatz**2)
        ehat1y = bhatx/np.sqrt(1-bhatz**2)
        ehat1z = 0*bhaty
        ehat2x = -bhatz*ehat1y 
        ehat2y = bhatz*ehat1x 
        ehat2z = bhatx*ehat1y-bhaty*ehat1x
        a = Pexx*ehat1x**2+Peyy*ehat1y**2+Pezz*ehat1z**2 + 2*(ehat1x*(Pexy*ehat1y+Pexz*ehat1z)+ehat1z*ehat1y*Peyz);
        d = Pexx*ehat2x**2+Peyy*ehat2y**2+Pezz*ehat2z**2 + 2*(ehat2x*(Pexy*ehat2y+Pexz*ehat2z)+ehat2z*ehat2y*Peyz);
        b = Pexx*ehat1x*ehat2x+Peyy*ehat1y*ehat2y+Pezz*ehat1z*ehat2z + Pexy*(ehat1x*ehat2y+ehat2x*ehat1y) \
            + Pexz*(ehat1x*ehat2z+ehat2x*ehat1z)+ Peyz*(ehat1z*ehat2y+ehat2z*ehat1y)
        Pperp1 = (a+d)/2+np.sqrt(((a-d)/2)**2+b**2) 
        Pperp2 = (a+d)/2-np.sqrt(((a-d)/2)**2+b**2) 
        Pp1p2 = 0*Ppar
        p1 = Pexx*ehat1x*bhatx+Peyy*ehat1y*bhaty+Pezz*ehat1z*bhatz + Pexy*(ehat1x*bhaty+bhatx*ehat1y) \
             + Pexz*(ehat1x*bhatz+bhatx*ehat1z)+ Peyz*(ehat1z*bhaty+bhatz*ehat1y);
        p2 = Pexx*ehat2x*bhatx+Peyy*ehat2y*bhaty+Pezz*ehat2z*bhatz + Pexy*(ehat2x*bhaty+bhatx*ehat2y) \
             + Pexz*(ehat2x*bhatz+bhatx*ehat2z)+ Peyz*(ehat2z*bhaty+bhatz*ehat2y);
        theta = .5*np.arcsin(2*b/(a+d))
        theta[np.isnan(theta)] = 0
        ct = np.cos(theta)
        st = np.sin(theta)

        fhat1x = ehat1x*ct - ehat2x*st 
        fhat1y = ehat1y*ct - ehat2y*st
        fhat1z = ehat1z*ct - ehat2z*st
        fhat2x = ehat2x*ct + ehat1x*st
        fhat2y = ehat2y*ct + ehat1y*st
        fhat2z = ehat2z*ct + ehat1z*st
        
        neUUep1p2 = neUUexx*fhat1x*fhat2x+neUUeyy*fhat1y*fhat2y+neUUezz*fhat1z*fhat2z + neUUexy*(fhat1x*fhat2y+fhat2x*fhat1y) \
                    + neUUexz*(fhat1x*fhat2z+fhat2x*fhat1z)+ neUUeyz*(fhat1z*fhat2y+fhat2z*fhat1y)
        neUUeparp1 = neUUexx*fhat1x*bhatx+neUUeyy*fhat1y*bhaty+neUUezz*fhat1z*bhatz + neUUexy*(fhat1x*bhaty+bhatx*fhat1y) \
                     + neUUexz*(fhat1x*bhatz+bhatx*fhat1z)+ neUUeyz*(fhat1z*bhaty+bhatz*fhat1y)
        neUUeparp2 = neUUexx*bhatx*fhat2x+neUUeyy*bhaty*fhat2y+neUUezz*bhatz*fhat2z + neUUexy*(bhatx*fhat2y+fhat2x*bhaty) \
                     + neUUexz*(bhatx*fhat2z+fhat2x*bhatz)+ neUUeyz*(bhatz*fhat2y+fhat2z*bhaty)
        neUUeperp1 = neUUexx*fhat1x**2+neUUeyy*fhat1y**2+neUUezz*fhat1z**2 + 2*(fhat1x*(neUUexy*fhat1y+neUUexz*fhat1z)+fhat1z*fhat1y*neUUeyz)
        neUUeperp2 = neUUexx*fhat2x**2+neUUeyy*fhat2y**2+neUUezz*fhat2z**2 + 2*(fhat2x*(neUUexy*fhat2y+neUUexz*fhat2z)+fhat2z*fhat2y*neUUeyz)
        
        Pparp1 = p1*np.cos(theta)-p2*np.sin(theta)
        Pparp2 = p2*np.cos(theta)+p1*np.sin(theta)
    
    
        Tpar = Ppar/ne
        Tperp = (Pperp1+Pperp2)/2/ne
        
        sz = np.shape(bx)
        Ibx = np.cumsum(bz[0,:])-(bz[0,:]+bz[0,0])/2
        Ibz = np.cumsum(bx, axis=0) - (np.matmul(np.ones((sz[0],1)),np.reshape(bx[0,:],(1,-1))) + bx)/2
        Psi = 2*(np.matmul(np.ones((sz[0],1)), np.reshape(Ibx,(1,-1))) - Ibz)*dx_de

        #Psi = (-np.cumsum(bx,axis=0)+np.cumsum(bz,axis=1))*dx_de;
        mfcont.update({'Ppar_ave':Ppar,'Pperp1_ave':Pperp1,'Pperp2_ave':Pperp2,'Pparp1_ave':Pparp1,
                       'Pparp2_ave':Pparp2,'Pp1p2_ave':Pp1p2,'Psi_ave':Psi,'Tpar_ave':Tpar,'Tperp_ave':Tperp,
                       'neUUepar_ave':neUUepar,'neUUeperp1_ave':neUUeperp1, 'neUUeperp2_ave':neUUeperp2, 'neUUeparp1_ave':neUUeparp1,
                       'neUUeparp2_ave':neUUeparp2,'neUUep1p2_ave':neUUep1p2})

        sio.savemat(matdir+ '/Slice_'+ str(slicenum)+ '.mat',mfcont)
    else:
        print('Slice ' + str(slicenum) + ' is missing averaged quantities!\n')

def processSliceavg(matdir, Slicenums):
    with open(matdir + '/../../info','r') as fp:
        info = fp.read()

    mime = int(float(info.split('mi/me =')[1].split()[0]))
    dx_de = float(info.split('dx/de =')[1].split()[0])
    wpewce = float(info.split('wpe/wce =')[1].split()[0])
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    dt = int(1/dtwci)

    #dt = mime*wpewce

    #for slicenum in Slicenums:
     #   process1slice(matdir,slicenum,dt,dx_de)
    slices = []
    for Slice in Slicenums:
        if os.path.isfile(matdir+'/../../data-ave/bx_'+str(int(dt*Slice)) +'.gda'):
            slices.append(Slice)
    pool = multiprocessing.Pool()
    pool.starmap(process1sliceavg,zip(repeat(matdir),slices,repeat(dt),repeat(dx_de)))
    FIXPSIAVG(matdir,slices)
    return 0

def FIXPSIAVG(matdir,slices):
    pool = multiprocessing.Pool()
    pool.starmap(fixPsi_ave,zip(repeat(matdir),slices))

def bundleslicesavg(basedir,slicenums):
    indir = basedir + '/Slices/'
    gdadir = basedir + '/data/'
    avdir = basedir + '/data-ave/'
    sliceavdir = indir +'SliceAve/'

    if not os.path.isdir(indir+'FullSlice/'):
        mkdir(indir + 'FullSlice/')
    
    avvars = []
    avts = []
    varlist = []
    avvarz = [f for f in listdir(sliceavdir) if  ((os.path.isfile(sliceavdir+f)) and (f!='info'))]
    avvarz = np.sort(avvarz)
    avlength = len(avvarz)
    varz = [f for f in listdir(indir) if ((os.path.isfile(indir+f)) and (f!='info'))]
    varz = np.sort(varz)
    length = len(varz)
    ts = []
    i = 0
    flag = True
    while flag:
        temp = varz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(ts))==0:
            ts.append(tnew)
        else:
            flag = False
    ts = np.sort(ts).astype(int)
    i = 0
    flag = True
    while flag:
        temp = avvarz[i].split('_ave_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(avts))==0:
            avts.append(tnew)
        else:
            flag = False
    avts = np.sort(avts).astype(int)
    j = -1
    for i in range(0,int(avlength/len(avts))):
        temp = avvarz[i*len(avts)].split('_')
        avvars.append(temp[0]+'_ave')

    tsdict = {}
    for i in range(0,(len(ts))):
        tsdict.update({int(ts[i]):i})
    avtsind = []
    for i in range(0,len(avts)-1):
        avtsind.append(tsdict[int(avts[i])])
    avtsind = np.array(avtsind) 
    avtsind = avtsind[np.where(avtsind>=np.min(slicenums))]
    avtsind = avtsind[np.where(avtsind<=np.max(slicenums))]
    #for slicenum in slicenums:
     #   bundler(indir,varlist,slicenum)
    pool = multiprocessing.Pool()
    pool.starmap(bundleravg,zip(repeat(indir),repeat(avvars),avtsind))
    return 0

def bundleravg(indir,varlist,slicenum):
    avmatdir = indir + '/SliceAve/'
    savedict = {}
    for i in range(0,len(varlist)):
        savedict.update(sio.loadmat(avmatdir+ '/'+ varlist[i]+ '_'+ str(slicenum)+ '.mat'))
        #savedict.update(sio.loadmat(indir+ '/' + varlist[i]+ str(slicenum)+ '.mat'))
    outfile = indir + 'FullSlice/Slice_'+ str(slicenum) + '.mat'
    dic = sio.loadmat(outfile)
    dic.update(savedict)
    sio.savemat(outfile, dic)

def populateslicesavg(varlist,ts,avts,nx,nz,avdir,savedir):
    tsdict = {}
    for i in range(0,(len(ts))):
        tsdict.update({int(ts[i]):i})
    avtsind = []
    for i in range(0,len(avts)):
        avtsind.append(tsdict[int(avts[i])])
    pool = multiprocessing.Pool()
    for var in varlist:
        pool.starmap(saveavgslice, zip(repeat(var),avts,avtsind,repeat(nx),repeat(nz),repeat(avdir),repeat(savedir)))
        #for slicenumM in range(0,len(ts)):
         #   saveslice(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir)
        #Parallel(n_jobs=-1)(delayed(saveslice)(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir) for slicenumM in range(0,len(ts)-1)) #now parallelized over time slice
    return 0

#def populateslicesvarpar(varlist,ts,nx,nz,gdadir,savedir):
 #   for slicenumM in range(0,len(ts)-1):
  #      Parallel(n_jobs=-1)(delayed(saveslice)(var,ts[slicenumM],slicenumM,nx,nz,gdadir,savedir) for var in varlist) #now parallelized over variable for smaller number of time slices
   # return 0

def makematsav(basedir):
    savedir = basedir + '/Slices/SliceAve/'
    gdadir = basedir + '/data/'
    avdir = basedir + '/data-ave/'

    if not os.path.isdir(savedir):
        mkdir(savedir)
    
    varlist = []
    varz = [f for f in listdir(gdadir) if ((os.path.isfile(gdadir+f)) and (f!='info'))]
    varz = np.sort(varz)
    length = len(varz)
    ts = []
    i = 0
    flag = True
    while flag:
        temp = varz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(ts))==0:
           ts.append(tnew)
        else:
            flag = False
    ts = np.sort(ts).astype(int)
    j = -1
    for i in range(0,int(length/len(ts))):
        temp = varz[i*len(ts)].split('_')
        varlist.append(temp[0])
    with open(basedir + '/info','r') as fp:
        info = fp.read()

    avvarlist = []
    avarz = [f for f in listdir(avdir) if ((os.path.isfile(avdir+f)) and (f!='info'))]
    avarz = np.sort(avarz)
    avlength = len(avarz)
    avts = []
    i = 0
    flag = True
    while flag:
        temp = avarz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(avts))==0:
           avts.append(tnew)
        else:
            flag = False
    avts = np.sort(avts).astype(int)
    j = -1
    for i in range(0,int(avlength/len(avts))):
        temp = avarz[i*len(avts)].split('_')
        avvarlist.append(temp[0])
    #print(avts)

    nsteps = int(info.split('num_step =')[1].split()[0])
    nx = int(float(info.split('nx =')[1].split()[0]))
    nz = int(float(info.split('nz =')[1].split()[0]))
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    nslices = int(nsteps*dtwci)
    fp.close()
    populateslicesavg(avvarlist,ts,avts,nx,nz,avdir,savedir)
    bundleslicesavg(basedir,range(0,nslices))
    processSliceavg(savedir+'../FullSlice/',range(0,nslices))

def contourvaluesav(indir, slicenums, xv, zv, Psivals, mime):
    ns = 1
    for slicenum in slicenums:
        if os.path.isfile(indir + '/../SliceAve/bx_ave_'+str(slicenum)+'.mat'):
            fdict = sio.loadmat(indir + '/Slice_' + str(slicenum) + '.mat')
            absB = np.sqrt(fdict['bx_ave']**2+fdict['by_ave']**2+fdict['bz_ave']**2)
            bhatx = fdict['bx_ave']/absB
            bhaty = fdict['by_ave']/absB
            bhatz = fdict['bz_ave']/absB

            ex = sp.ndimage.filters.gaussian_filter(fdict['ex_ave'],ns,mode='constant')
            ey = sp.ndimage.filters.gaussian_filter(fdict['ey_ave'],ns,mode='constant')
            ez = sp.ndimage.filters.gaussian_filter(fdict['ez_ave'],ns,mode='constant')

            Epar = ex*bhatx+ey*bhaty+ez*bhatz

            splineB = sp.interpolate.RectBivariateSpline(xv,zv,absB.T)
            splineE = sp.interpolate.RectBivariateSpline(xv,zv,Epar.T)
            splineBg = sp.interpolate.RectBivariateSpline(xv,zv,bhaty.T)
            splinebx = sp.interpolate.RectBivariateSpline(xv,zv,bhatx.T)
            splinebz = sp.interpolate.RectBivariateSpline(xv,zv,bhatz.T)
            newdict = {}
            for j in range(0,len(Psivals)): 
                temp = plt.contour(xv,zv,fdict['Psi_ave'],Psivals).collections[j].get_paths()[0].vertices
                x = temp[:,0]
                z = temp[:,1]
                #print(xv,zv,x,z)
                B = splineB.ev(x,z)
                B[B!=B] = 1
                    
                epar = splineE.ev(x,z)
                epar[epar!=epar]  = 0
                dx = (np.append(0,np.diff(x))+np.append(np.diff(x),0))/2
                dz = (np.append(0,np.diff(z))+np.append(np.diff(z),0))/2

                bg = splineBg.ev(x,z)
                bg[bg!=bg] = 0
                bx = splinebx.ev(x[0],z[0])
                bz = splinebz.ev(x[0],z[0])
                parsign = np.sign(bx*dx[0]+bz*dz[0])
                dl =  parsign*np.sqrt((dx**2+dz**2)*(1/(1-bg**2)))*np.sqrt(mime)
                l = np.cumsum(dl)
                phipar = -np.cumsum(dl*epar)
                midind = np.argmin(np.abs(z))
                l = l - l[midind]
                phipar = phipar - phipar[midind]
                newdict.update({'x'+str(j):x,'z'+str(j):z,'l'+str(j):l,'B'+str(j):B,'phipar'+str(j):phipar})
                dic = sio.loadmat(indir + '/Slice_'+ str(slicenum)+ '.mat')
                dic.update({'xv':xv,'zv':zv,'Psivals':Psivals})
                dic.update(newdict)
                sio.savemat(indir + '/Slice_'+ str(slicenum)+ '.mat', dic)
    
def fulldJdt(indir,t0,t1,nPsi,vA,wpewce,mime):
    m = 1
    e = -1
    dtwpe = wpewce*mime
    slicevals = sio.loadmat(indir+'/Slice_'+str(t0)+'.mat')
    Jmat0 = slicevals['Jmat']
    mumat0 = slicevals['mumat']
    slicevals = sio.loadmat(indir+'/Slice_'+str(t1)+'.mat')
    Jmat1 = slicevals['Jmat']
    mumat1 = slicevals['mumat']
    nU,nmu,_ = np.shape(Jmat0)
    mumax = np.zeros((nU,1))
    mumin = np.zeros((nU,1))
    J0 = np.zeros(np.shape(Jmat0[:,:,0]))
    J1 = np.zeros(np.shape(Jmat1[:,:,0]))
    dJdt = np.zeros(np.shape(Jmat0))
    mumat = np.zeros(np.shape(Jmat0))
    for j in range(0,np.size(nPsi)):
        for i in range(0,nU):
            k = nPsi[j]
            mumax = min((mumat0[i,nmu-1,k],mumat1[i,nmu-1,k]))
            mumin = max((mumat0[i,0,k],mumat1[i,0,k]))
            mu = np.linspace(mumin,mumax,nmu)
            mumat[i,:,j] = mu;
            J0[i,:] = np.interp(mu,mumat0[i,:,k],Jmat0[i,:,k])
            J1[i,:] = np.interp(mu,mumat1[i,:,k],Jmat1[i,:,k])
        dJdt[:,:,j] = (J1-J0)/((t1-t0)*dtwpe)
    sio.savemat(indir + '../dJdt_'+str(t0)+'.mat', {'dJdt':dJdt,'mumat_d':mumat})

def makePsiGIF(fullslicedir,slicenums,contourlevels,xv,zv):
    fig = plt.figure()
    psigif = FuncAnimation(fig, plotPsiCont, slicenums , fargs = (fullslicedir,contourlevels,xv,zv),interval=50,repeat=False)
    psigif.save(fullslicedir+'../../PsiGIF.gif', dpi=75, writer='imagemagick')
    
def plotPsiCont(slicenum,fullslicedir,contourlevels,xv,zv):
    a = sio.loadmat(fullslicedir+'/Slice_'+str(slicenum)+'.mat')
    Psi = a['Psi']
    #print(np.shape(xv),np.shape(zv),np.shape(Psi))
    ax = plt.contour(xv,zv,Psi,levels=contourlevels)
    plt.xlabel('x (d_i)')
    plt.ylabel('z (d_i)')
    plt.title('Psi: Slice '+ str(slicenum))
    plt.draw()
    print('Finished Slice '+str(slicenum)+'.')
    return ax

def populatePsi(fullslicedir,slicenums,contourlevels,xv,zv):
    if not os.path.isdir(fullslicedir+'../../Images/'):
        mkdir(fullslicedir+'../../Images')
    pool = multiprocessing.Pool()
    pool.starmap(makePsiSlice, zip(slicenums,repeat(fullslicedir),repeat(contourlevels),repeat(xv),repeat(zv)))
    makeGIF(fullslicedir+'../../Images/','Psi', slicenums, '.png')

def makePsiSlice(slicenum,fullslicedir,contourlevels,xv,zv):
    fig = plt.figure()
    a = sio.loadmat(fullslicedir+'/Slice_'+str(slicenum)+'.mat')
    Psi = a['Psi']
    ax = plt.contour(xv,zv,Psi,levels=contourlevels)
    plt.xlabel('x (d_i)')
    plt.ylabel('z (d_i)')
    plt.title('Psi: Slice '+ str(slicenum))
    plt.draw()
    plt.savefig(fullslicedir +'../../Images/Psi_'+str(slicenum)+'.png')

def makeGIF(imdir, basename, slicenums, imageext):
    images = [(imdir + basename + '_' + str(index) + imageext) for index in slicenums]
    filename = imdir+'../'+basename+'.gif'
    with open(imdir + basename+'_list.txt','w') as fil:
        for item in images:
            fil.write("%s\n" % item)
    os.chdir(imdir)
    os.system('convert @'+basename+'_list.txt '+filename)
    #clip = mpy.ImageSequenceClip(images,fps=10)
    #clip.write_gif(filename,fps=10,program='imagemagick')
    #writeGif(filename,images,duration=0.1)
    #with imageio.get_writer(imdir+'../'+basename+'.gif', mode='I') as writer:
     #   for index in range(startindex,endindex):
      #      filename = basename + '_' + str(index) + imageext
       #     image = imageio.imread(filename)
        #    writer.append_data(image)
def populatePcolor(plotvar,fullslicedir,slicenums,xv,zv):
    if not os.path.isdir(fullslicedir+'../../Images/'):
        mkdir(fullslicedir+'../../Images/')
    pool = multiprocessing.Pool()
    pool.starmap(makePcolorSlice, zip(repeat(plotvar),slicenums,repeat(fullslicedir),repeat(xv),repeat(zv)))
    makeGIF(fullslicedir+'../../Images/',plotvar, slicenums, '.jpeg')

def makePcolorSlice(plotvar,slicenum,fullslicedir,xv,zv):
    fig = plt.figure()
    a = sio.loadmat(fullslicedir+'/Slice_'+str(slicenum)+'.mat')
    A = a[plotvar]
    ax = plt.pcolormesh(xv,zv,A)
    plt.colorbar()
    plt.xlabel('x (d_i)')
    plt.ylabel('z (d_i)')
    plt.title(plotvar + ': Slice '+ str(slicenum))
    plt.draw()
    plt.savefig(fullslicedir +'../../Images/'+plotvar+'_'+str(slicenum)+'.jpeg')    

def populatePcolorPsi(plotvar,fullslicedir,slicenums,contourlevels,xv,zv):
    if not os.path.isdir(fullslicedir+'../../Images/'):
        mkdir(fullslicedir+'../../Images/')
    pool = multiprocessing.Pool()
    pool.starmap(makePcolorPsiSlice, zip(repeat(plotvar),slicenums,repeat(contourlevels),repeat(fullslicedir),repeat(xv),repeat(zv)))
    makeGIF(fullslicedir+'../../Images/',plotvar+'_cont', slicenums, '.jpeg')

def makePcolorPsiSlice(plotvar,slicenum,contourlevels,fullslicedir,xv,zv):
    fig = plt.figure()
    a = sio.loadmat(fullslicedir+'/Slice_'+str(slicenum)+'.mat')
    A = a[plotvar]
    Psi = a['Psi']
    ax = plt.pcolormesh(xv,zv,A)
    plt.colorbar()
    plt.hold(True)
    plt.contour(xv,zv,Psi,levels=contourlevels,colors='black',linewidths=2)
    plt.xlabel('x (d_i)')
    plt.ylabel('z (d_i)')
    plt.title(plotvar + ': Slice '+ str(slicenum))
    plt.draw()
    plt.savefig(fullslicedir +'../../Images/'+plotvar+'_cont_'+str(slicenum)+'.jpeg')

def makePsivals(slicedir,slicenum,xv,zv, xmin, xmax, nPsi):
    xs = np.linspace(xmin,xmax,nPsi)
    Psi = sio.loadmat(slicedir + '/Slice_' + str(slicenum) + '.mat')['Psi']
    splinePsi = sp.interpolate.RectBivariateSpline(xv,zv,Psi.T)
    Psivals = splinePsi.ev(xs,np.zeros(np.shape(xs)))
    return np.sort(Psivals)

def contourvaluesavx(indir, slicenums, xmin, xmax, xv, zv, nPsi, mime):
    pool = multiprocessing.Pool()
    pool.starmap(contourvaluesavx1slice, zip(repeat(indir),slicenums,repeat(xmin),repeat(xmax),repeat(xv),repeat(zv),repeat(nPsi),repeat(mime))) 

def contourvaluesavx1slice(indir,slicenum,xmin,xmax,xv,zv,nPsi,mime):
    ns = 1
    if not (slicenum % 2):
        Psivals = makePsivals(indir,slicenum,xv,zv,xmin,xmax,nPsi)
    else:
        Psivals = makePsivals(indir,slicenum-1,xv,zv,xmin,xmax,nPsi)
    if os.path.isfile(indir + '/../SliceAve/bx_ave_'+str(slicenum)+'.mat'):
        fdict = sio.loadmat(indir + '/Slice_' + str(slicenum) + '.mat')
        absB = np.sqrt(fdict['bx_ave']**2+fdict['by_ave']**2+fdict['bz_ave']**2)
        bhatx = fdict['bx_ave']/absB
        bhaty = fdict['by_ave']/absB
        bhatz = fdict['bz_ave']/absB

        ex = sp.ndimage.filters.gaussian_filter(fdict['ex_ave'],ns,mode='constant')
        ey = sp.ndimage.filters.gaussian_filter(fdict['ey_ave'],ns,mode='constant')
        ez = sp.ndimage.filters.gaussian_filter(fdict['ez_ave'],ns,mode='constant')

        Epar = ex*bhatx+ey*bhaty+ez*bhatz

        splineB = sp.interpolate.RectBivariateSpline(xv,zv,absB.T)
        splineE = sp.interpolate.RectBivariateSpline(xv,zv,Epar.T)
        splineBg = sp.interpolate.RectBivariateSpline(xv,zv,bhaty.T)
        splinebx = sp.interpolate.RectBivariateSpline(xv,zv,bhatx.T)
        splinebz = sp.interpolate.RectBivariateSpline(xv,zv,bhatz.T)
        newdict = {}
        for j in range(0,nPsi): 
            temp = plt.contour(xv,zv,fdict['Psi_ave'],Psivals).collections[j].get_paths()[0].vertices
            x = temp[:,0]
            z = temp[:,1]
            #print(xv,zv,x,z)
            B = splineB.ev(x,z)
            B[B!=B] = 1
            
            epar = splineE.ev(x,z)
            epar[epar!=epar]  = 0
            dx = (np.append(0,np.diff(x))+np.append(np.diff(x),0))/2
            dz = (np.append(0,np.diff(z))+np.append(np.diff(z),0))/2

            bg = splineBg.ev(x,z)
            bg[bg!=bg] = 0
            bx = splinebx.ev(x[0],z[0])
            bz = splinebz.ev(x[0],z[0])
            parsign = np.sign(bx*dx[0]+bz*dz[0])
            dl =  parsign*np.sqrt((dx**2+dz**2)*(1/(1-bg**2)))*np.sqrt(mime)
            l = np.cumsum(dl)
            phipar = -np.cumsum(dl*epar)
            midind = np.argmin(np.abs(z))
            l = l - l[midind]
            phipar = phipar - phipar[midind]
            newdict.update({'x'+str(j):x,'z'+str(j):z,'l'+str(j):l,'B'+str(j):B,'phipar'+str(j):phipar})
            dic = sio.loadmat(indir + '/Slice_'+ str(slicenum)+ '.mat')
            dic.update({'xv':xv,'zv':zv,'Psivals':Psivals})
            dic.update(newdict)
            sio.savemat(indir + '/Slice_'+ str(slicenum)+ '.mat', dic)

def contourvaluesx(indir, slicenums, xmin, xmax, xv, zv, nPsi, mime):
    pool = multiprocessing.Pool()
    pool.starmap(contourvaluesx1slice, zip(repeat(indir),slicenums,repeat(xmin),repeat(xmax),repeat(xv),repeat(zv),repeat(nPsi),repeat(mime))) 

def contourvaluesx1slice(indir,slicenum,xmin,xmax,xv,zv,nPsi,mime):
    ns = 10
    if not (slicenum % 2):
        Psivals = makePsivals(indir,slicenum,xv,zv,xmin,xmax,nPsi)
    else:
        Psivals = makePsivals(indir,slicenum-1,xv,zv,xmin,xmax,nPsi)
    fdict = sio.loadmat(indir + '/Slice_' + str(slicenum) + '.mat')
    absB = np.sqrt(fdict['bx']**2+fdict['by']**2+fdict['bz']**2)
    bhatx = fdict['bx']/absB
    bhaty = fdict['by']/absB
    bhatz = fdict['bz']/absB

    ex = sp.ndimage.filters.gaussian_filter(fdict['ex'],ns,mode='constant')
    ey = sp.ndimage.filters.gaussian_filter(fdict['ey'],ns,mode='constant')
    ez = sp.ndimage.filters.gaussian_filter(fdict['ez'],ns,mode='constant')

    Epar = ex*bhatx+ey*bhaty+ez*bhatz

    splineB = sp.interpolate.RectBivariateSpline(xv,zv,absB.T)
    splineE = sp.interpolate.RectBivariateSpline(xv,zv,Epar.T)
    splineBg = sp.interpolate.RectBivariateSpline(xv,zv,bhaty.T)
    splinebx = sp.interpolate.RectBivariateSpline(xv,zv,bhatx.T)
    splinebz = sp.interpolate.RectBivariateSpline(xv,zv,bhatz.T)
    newdict = {}
    for j in range(0,nPsi): 
        temp = plt.contour(xv,zv,fdict['Psi'],Psivals).collections[j].get_paths()[0].vertices
        x = temp[:,0]
        z = temp[:,1]
        B = splineB.ev(x,z)
        B[B!=B] = 1
            
        epar = splineE.ev(x,z)
        epar[epar!=epar]  = 0
        dx = (np.append(0,np.diff(x))+np.append(np.diff(x),0))/2
        dz = (np.append(0,np.diff(z))+np.append(np.diff(z),0))/2

        bg = splineBg.ev(x,z)
        bg[bg!=bg] = 0
        bx = splinebx.ev(x[0],z[0])
        bz = splinebz.ev(x[0],z[0])
        parsign = np.sign(bx*dx[0]+bz*dz[0])
        dl =  parsign*np.sqrt((dx**2+dz**2)*(1/(1-bg**2)))*np.sqrt(mime)
        l = np.cumsum(dl)
        phipar = -np.cumsum(dl*epar)
        midind = np.argmin(np.abs(z))
        l = l - l[midind]
        phipar = phipar - phipar[midind]
        newdict.update({'x'+str(j):x,'z'+str(j):z,'l'+str(j):l,'B'+str(j):B,'phipar'+str(j):phipar})
        dic = sio.loadmat(indir + '/Slice_'+ str(slicenum)+ '.mat')
        dic.update({'xv':xv,'zv':zv,'Psivals':Psivals})
        dic.update(newdict)
        sio.savemat(indir + '/Slice_'+ str(slicenum)+ '.mat', dic)

def makematsserial(basedir):
    savedir = basedir + '/Slices/'
    gdadir = basedir + '/data/'

    if not os.path.isdir(savedir):
        mkdir(savedir)
    
    varlist = []
    varz = [f for f in listdir(gdadir) if ((os.path.isfile(gdadir+f)) and (f!='info'))]
    varz = np.sort(varz)
    length = len(varz)
    ts = []
    i = 0
    flag = True
    while flag:
        temp = varz[i].split('_')
        i = i + 1
        temp2 = temp[1].split('.')
        tnew = int(temp2[0])
        if np.sum(tnew==np.array(ts))==0:
           ts.append(tnew)
        else:
            flag = False
    ts = np.sort(ts).astype(int)
    j = -1
    for i in range(0,int(length/len(ts))):
        temp = varz[i*len(ts)].split('_')
        varlist.append(temp[0])
    with open(basedir + '/info','r') as fp:
        info = fp.read()

    nsteps = int(info.split('num_step =')[1].split()[0])
    nx = int(float(info.split('nx =')[1].split()[0]))
    nz = int(float(info.split('nz =')[1].split()[0]))
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    nslices = int(nsteps*dtwci);
   # nx = 2560;
   # nz = 2560;
    fp.close()
    populateslices(varlist,ts,nx,nz,gdadir,savedir)
    bundleslices(basedir,range(0,len(ts)))
    processSliceSerial(savedir+'/FullSlice/',range(0,len(ts)))

def processSliceSerial(matdir, Slicenums):
    with open(matdir + '/../../info','r') as fp:
        info = fp.read()

    mime = int(float(info.split('mi/me =')[1].split()[0]))
    dx_de = float(info.split('dx/de =')[1].split()[0])
    wpewce = float(info.split('wpe/wce =')[1].split()[0])
    dtwci = float(info.split('dt*wci =')[1].split()[0])
    #v_A = float(info.split('v_A')[1].split('=')[1].split()[0])
    dt0 = wpewce*mime
    dt = int(1/dtwci)
    Edrive = float(info.split('edrive =')[1].split()[0])
    taudrive = float(info.split('tdrive =')[1].split()[0])

    #taudrive = dt0*tauwci
    #Edrive = v_A/wpewce*edrivefactor
    #dt = mime*wpewce

    for slicenum in Slicenums:
        process1slice(matdir,slicenum,dt,dx_de)
    fixPsi(matdir,Slicenums,dt0, Edrive, taudrive)
    return 0

def contourvaluesxserial(indir, slicenums, xmin, xmax, xv, zv, nPsi, mime):
    for slicenum in slicenums:
        contourvaluesx1slice(indir,slicenum,xmin,xmax,xv,zv,nPsi,mime)

def loadsliceserial(q,slicenum,nx,nz,datadir):
    q = ''.join(q.split('_ave'))
    casedict = {
        'Pexx':'pe-xx',
        'Pexy':'pe-xy',
        'Pexz':'pe-xz',
        'Peyy':'pe-yy',
        'Peyz':'pe-yz',
        'Pezz':'pe-zz',
        'Pixx':'pi-xx',
        'Pixy':'pi-xy',
        'Pixz':'pi-xz',
        'Piyy':'pi-yy',
        'Piyz':'pi-yz',
        'Pizz':'pi-zz',
        'neUUexx':'neuue-xx',
        'neUUexy':'neuue-xy',
        'neUUexz':'neuue-xz',
        'neUUeyy':'neuue-yy',
        'neUUeyz':'neuue-yz',
        'neUUezz':'neuue-zz'
    }
    q = casedict.get(q,q)
    
    with open(datadir+'/'+ q + '.gda','rb') as fd:
        for i in range(0,slicenum+1):
            a = np.reshape(np.fromfile(fd,dtype='f4',count = nx*nz),(nz,nx))
            np.fromfile(fd,dtype='f4',count=2)
        return a

def savesliceserial(var,slicenum,nx,nz,gdadir,savedir):
    casedict = {
        'pe-xx':'Pexx',
        'pe-xy':'Pexy',
        'pe-xz':'Pexz',
        'pe-yy':'Peyy',
        'pe-yz':'Peyz',
        'pe-zz':'Pezz',
        'pi-xx':'Pixx',
        'pi-xy':'Pixy',
        'pi-xz':'Pixz',
        'pi-yy':'Piyy',
        'pi-yz':'Piyz',
        'pi-zz':'Pizz'
    }
    var = casedict.get(var,var)
    outfile = savedir + '/' + var +'_' + str(slicenum)  +'.mat'
    #outfile  = savedir + '/' + var + str(slicenumM) + '.mat'

    if not (os.path.isfile(outfile)):
        out = loadsliceserial(var,slicenum,nx,nz,gdadir)
        out = out[::2,:]+out[1:out.shape[0]:2,:] #decimation steps
        out = out[:,::2]+out[:,1:out.shape[1]:2]
        sio.savemat(outfile,{var:out/4})
    return 0
