from matfunctions import *
import numpy as np
#basedir = '/net/scratch3/wetherton/bg1.75-mime100/'
#basedir = '/net/scratch3/wetherton/testavg4/'
basedir = '/net/scratch3/wetherton/dJdt/Run4/'
#basedir = '/lustre/scratch3/turquoise/daughton/Blake-time-average/'
makemats(basedir)
print('Finished makemats')
makematsav(basedir)
print('Finished makematsav')
#contourvaluesavx(basedir + 'Slices/FullSlice/',range(0,300), 10, 40,np.linspace(0,100,3150),np.linspace(-10,10,630), 20, 100)
#contourvaluesav(basedir+'Slices/FullSlice/',(220,221),np.linspace(0,100,3150),np.linspace(-10,10,630),np.linspace(220,250,10),100)
#print('Finished Contour')
#calcJ(basedir+'Slices/FullSlice/',(20,21),1.5,100,50)
#contourvaluesav(basedir+'Slices/FullSlice/',(10,11),np.linspace(0,40,2520),np.linspace(-10,10,1260),[9.5, 10, 10.5, 11.0, 11.5],400)
#calcJ(basedir+'Slices/FullSlice/',(220,221),0.5,100,100)
#print('Calculated J')
#alpha1,alpha2 =calcalpha12(basedir+'Slices/FullSlice/',220,221,range(0,20),.07071068,2.0,100)
#print(alpha1,alpha2)
#fulldJdt(basedir+'Slices/FullSlice/',220,221,range(0,20),.07071068,2.0,100)


