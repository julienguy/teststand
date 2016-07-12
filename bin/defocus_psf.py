#!/usr/bin/env python

import os,sys
import astropy.io.fits as pyfits
import numpy as np
import pylab
import scipy.interpolate
if not "DESIMODEL" in os.environ :
    print "need DESIMODEL env. variable"
    sys.exit(12)

psffile="%s/data/specpsf/psf-r.fits"%(os.environ["DESIMODEL"])
desimodel_psf=pyfits.open(psffile)
desimodel_psf.info()
head=desimodel_psf[0].header
print head
desimodel_CCDPIXSZ=head["CCDPIXSZ"]*1e3 # conversion from mm to um
desimodel_PIXSIZE=head["PIXSIZE"]*1e3 # conversion from mm to um

print "desimodel_CCDPIXSZ",desimodel_CCDPIXSZ,"um"
print "desimodel_PIXSIZE",desimodel_PIXSIZE,"um"

defocus_PIXSIZE=0

defocus_psf_pos=[]
defocus_psf_wave=[]
defocus_psf_filename=[]
defocus_psf_3dimage=[]
for f in os.listdir(".") :
    if f.find('DESI_FOCUSTEST_RED_')==0 :
        pos=float(f.split("fx=")[1].split("mm")[0]) 
        wave=float(f.split("_wl=")[1].split("um")[0])*1e4 # A
        print f,pos,wave
        defocus_psf_pos.append(pos)
        defocus_psf_wave.append(wave)
        defocus_psf_filename.append(f)
        h=pyfits.open(f)
        
        if defocus_PIXSIZE==0 :
            defocus_PIXSIZE=h[0].header["PIXEL"]
        else :
            if defocus_PIXSIZE != h[0].header["PIXEL"] :
                print "warning DIFFERENT PIXEL SIZE"
                sys.exit(12)
        
        defocus_psf_3dimage.append(h[0].data.astype(float))

head=pyfits.open(defocus_psf_filename[0])[0].header
defocus=np.zeros((31))
for i in range(31) :
    defocus[i]=float(head["DEFOC%02d"%i])
print "defocus values (in mm) :",defocus
    
t_wave=np.unique(defocus_psf_wave)
t_pos=np.unique(defocus_psf_pos)
print "defocus waves = ",t_wave
print "defocus positions = ",t_pos
print "desimodel waves = ",desimodel_psf["SPOTWAVE"].data

spots=desimodel_psf["SPOTS"].data

for defocus_index in range(len(defocus)) :
    print "doing defocus values (in mm) :",defocus[defocus_index]

    for i in range(spots.shape[0]) :
        for j in range(spots.shape[1]) :
            spot_wave=desimodel_psf["SPOTWAVE"].data[j]
            spot_x_pix=desimodel_psf["SPOTX"].data[i,j]
            spot_y_pix=desimodel_psf["SPOTY"].data[i,j]
            
            
            # find nearest wave
            defocus_w = t_wave[np.argmin(np.abs(t_wave-spot_wave))]
            # find nearest pos
            spot_pos = (spot_x_pix-2050)*desimodel_CCDPIXSZ # mm
            defocus_pos = t_pos[np.argmin(np.abs(t_pos-spot_pos))]
            # find img            
            defoc_psf_index=np.where((defocus_psf_wave==defocus_w)&(defocus_psf_pos==defocus_pos))[0]
            
            print i,j,"w=",spot_wave,"pos=",spot_pos,"->","w=",defocus_w,"pos=",defocus_pos,"index=",defoc_psf_index

            plot=False

            if plot :
                fig=pylab.figure()
                pl0=pylab.subplot(2,2,1)
                pl1=pylab.subplot(2,2,2)
                pl2=pylab.subplot(2,2,3)
                pl3=pylab.subplot(2,2,4)
                
            spot=spots[i,j]
            dspot=defocus_psf_3dimage[defoc_psf_index][defocus_index]

            sx=np.arange(spot.shape[0])
            sy=np.arange(spot.shape[1])                
            sxcen=np.sum(sx*np.sum(spot,axis=1))/float(np.sum(spot))
            sycen=np.sum(sy*np.sum(spot,axis=0))/float(np.sum(spot))
            sx=(sx-sxcen)*desimodel_PIXSIZE/desimodel_CCDPIXSZ
            sy=(sy-sycen)*desimodel_PIXSIZE/desimodel_CCDPIXSZ


            dx=np.arange(dspot.shape[0])
            dy=np.arange(dspot.shape[1])             
            dxcen=np.sum(dx*np.sum(dspot,axis=1))/float(np.sum(dspot))
            dycen=np.sum(dy*np.sum(dspot,axis=0))/float(np.sum(dspot))
            dx=(dx-dxcen)*defocus_PIXSIZE/desimodel_CCDPIXSZ
            dy=(dy-dycen)*defocus_PIXSIZE/desimodel_CCDPIXSZ
                
                
            if plot :
                pl2.plot(dx,dspot[:,int(dycen)]/np.max(dspot))
                pl2.plot(sx,spot[:,int(sycen)]/np.max(spot))
                pl2.set_xlabel("x (CCD pix)")
                pl3.plot(dy,dspot[int(dxcen),:]/np.max(dspot))
                pl3.plot(sy,spot[int(sxcen),:]/np.max(spot))
                pl3.set_xlabel("y (CCD pix)")
                    
            # and now we try to interpolate all of this
            f = scipy.interpolate.interp2d(dy,dx, dspot, kind='linear')
            idspot = f(sy,sx)
                
            if plot :
                pl2.plot(sx,idspot[:,int(sycen)]/np.max(idspot))
                pl3.plot(sy,idspot[int(sxcen),:]/np.max(idspot))                
                pl0.imshow(idspot,extent=(sx[0],sx[-1],sy[0],sy[-1]),interpolation="nearest",origin=0)
                pl1.imshow(spot,extent=(sx[0],sx[-1],sy[0],sy[-1]),interpolation="nearest",origin=0)
                pylab.show()
                sys.exit(0)

            # replace the spot
            spots[i,j] = idspot

    # end of loop on i,j
    desimodel_psf["SPOTS"].data = spots
    ofilename="psf-r_defoc_%04.3fmm.fits"%defocus[defocus_index]
    desimodel_psf.writeto(ofilename,clobber=True)
    print "wrote",ofilename








