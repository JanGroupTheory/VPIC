from matfunctions_memory import *
import numpy as np
#basedir = '/net/scratch3/wetherton/bg1.75-mime100/'
#basedir = '/net/scratch3/wetherton/testavg4/'
basedir = '/cygdrive/g/dJdt/Run1/'
#basedir = '/lustre/scratch3/turquoise/daughton/Blake-time-average/'
#makemats(basedir)
#print('Finished makemats')
#makematsav(basedir)
#print('Finished makematsav')
#contourvaluesavx(basedir + 'Slices/SliceAve/',range(0,95), 85.1, 160, np.linspace(0,320,20480),np.linspace(-15,15,1920), 50, 400)
#contourvaluesav(basedir+'Slices/FullSlice/',(220,221),np.linspace(0,100,3150),np.linspace(-10,10,630),np.linspace(220,250,10),100)
#print('Finished Contour')
calcJ(basedir+'Slices/FullSlice/',(80,81),0.2,100,99)
#contourvaluesav(basedir+'Slices/FullSlice/',(10,11),np.linspace(0,40,2520),np.linspace(-10,10,1260),[9.5, 10, 10.5, 11.0, 11.5],400)
#calcJ(basedir+'Slices/FullSlice/',(220,221),0.5,100,100)
#print('Calculated J')
calcalpha12(basedir+'Slices/FullSlice/',80,81,range(0,50),.03535534,2.0,400)
#print(alpha1,alpha2)
fulldJdt(basedir+'Slices/FullSlice/',80,81,range(0,50),.03535534,2.0,400)


