#!/usr/bin/env python
##################################################################################
##################################################################################
# get_z_interactive.py By Nathaniel R. Ross, UCLA, nross@astro.ucla.edu
# usage: python get_z_interactive.py line_list
# Reads in list of emission lines from the WISP survey HST/WFC3/IR Grisms and 
# plots the spectra, iterates through the detected emission lines, allowing the 
# user to make line identifications or reject spurious lines quickly and 
# efficiently.
#
# Version 1.0 updates the program to look for wavelength-calibrated 2d grism 
# stamps first. Also, will now look for default line list name and make 
# linelistfile an optional parameter. Furthermore, I have added the option 
# to save the incomplete line list to file as you progress, but default is NOT to 
# do this.
# Optional arguments for go() are 1. linelistfile (String) - path to line list file
#                                 2. save_temp (Bool) - Save progress in 
#                                    linelist/Par???lines_with_redshifts.incomplete.dat
#                                 3. recover_temp (Bool) - Recover progress from previous 
#                                    session using the .incomplete.dat files
#
# ** Major change in version 1.0:
# We now use xpa instead of IRAF to display the 2d grism stamps and full 2d 
# direct images (instead of cutout stamp). Reason for this change was the desire
# to show the grism stamps that have a wavelength solution applied. When using 
# the IRAF display command, ds9 would fail to recognize this coordinate system.
# The XPA system can be downloaded here: http://hea-www.harvard.edu/RD/xpa/index.html
#
##################################################################################
import os
import distutils
import numpy as np
import scipy
import pylab as plt
from scipy.interpolate import spline
from distutils.sysconfig import *
import sys
import astropy.io.fits as pyfits

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def ret_z_lines (linelams,redshift,minlam,maxlam,minval,maxval):
    retx=[]
    rety=[]
    for j in range(len(linelams)):
        if (linelams[j]*(1+redshift))>=minlam and (linelams[j]*(1+redshift))<=maxlam:
            retx.append(linelams[j]*(1+redshift))
            retx.append(linelams[j]*(1+redshift))
            rety.append(minval)
            rety.append(maxval)
    return [retx,rety]

def getzeroorders (zeroorderpath,g='G141',magcut=23.5): # MB: changed from 23.0
    zop=open(zeroorderpath,'r')
    zox=[]
    zoy=[]
    zoid=[]
    zmag=[]
    for line in zop:
        if len(line)>60:
            linesplit = line.split()
            zox.append(float(linesplit[1][0:-1]))
            #zox.append(float(line[7:14]))   #### old syntax reg file 
            zoy.append(float(linesplit[2][0:-3])) 
            #zoid.append(int( linesplit[-2][-1]))  #### was only getting 1 digit
            zoid.append(int(linesplit[-2].split('{')[-1]))
            zmag.append(float(linesplit[-1][1:-2]))  ### get mag from reg file

            #zoy.append(float(line[16:23]))  ### old synatx reg file 
            #ll=len(line)
            #indend=ll-2
            #ind1=indend-6
            #for i in range(7):
            #    if line[ind1+i]!=' ' and line[ind1+i]!='{':
            #        indstart=ind1+i
            #        break
            #if indstart<indend:
            #    zoid.append(int(line[indstart:indend]))
            #else:
            #    zoid.append(int(line[indend]))
    zop.close()
    ### don't need the get the catalogs anymore ###
    ### getting magnitudes from region file ###
#    p_cat110='../DATA/DIRECT_GRISM/fin_F110.cat'
#    p_cat140='../DATA/DIRECT_GRISM/fin_F140.cat'
#    p_cat160='../DATA/DIRECT_GRISM/fin_F160.cat'
#    if (g=='G141' and os.path.exists(p_cat140)==0 and os.path.exists(p_cat160)==0):
#        print 'Cannot find H-band source catalog to determine magnitude cut on G141 Zero-order reg.'
#        return (zox,zoy,zoid)
#    elif g=='G141' and os.path.exists(p_cat140)==1:
#        print 'Using F140W source catalog to determine magnitude cut on G141 Zero-order reg.'
#        p=p_cat140
#    elif g=='G141' and os.path.exists(p_cat160)==1:
#        print 'Using F160W source catalog to determine magnitude cut on G141 Zero-order reg.'
#        p=p_cat160
#    if g!='G141' and os.path.exists(p_cat110)==0 and os.path.exists(p_cat140)==1:
#        print 'No F110W source catalog. Using F140W catalog for magnitude cut on G102 Zero-order reg.'
#        p=p_cat140
#    elif g!='G141' and os.path.exists(p_cat110)==0 and os.path.exists(p_cat160)==1:
#        print 'No F110W source catalog. Using F160W catalog for magnitude cut on G102 Zero-order reg.'
#        p=p_cat160
#    elif g!='G141' and os.path.exists(p_cat110)==1:
#        print 'Using F110W source catalog to determine magnitude cut on G102 Zero-order reg.'
#        p=p_cat110
#    elif g!='G141':
#        print 'No source catalog found to determine magnitude cut on G102 Zero-order reg.'
#        return (zox,zoy,zoid)
#    catdat=np.genfromtxt(p,dtype=np.str)
#    catid=np.array(catdat[0:,1],dtype=np.int)
#    catmags=np.array(catdat[0:,12],dtype=np.float)
#    filt=catmags<=magcut
    zoid=np.array(zoid)
    zoy=np.array(zoy)
    zox=np.array(zox)
    zmag=np.array(zmag)
    cond = (zmag <= magcut)
    return zox[cond],zoy[cond],zoid[cond]
#    zx,zy,zi=[],[],[]
#    for id_mag in catid[filt]:
#        zfilt=zoid==id_mag
#        if len(zoid[zfilt])==1:
#            zx.append(zox[zfilt])
#            zy.append(zoy[zfilt])
#            zi.append(zoid[zfilt])
#    return (zx,zy,zi)


def getfirstorders (firstorderpath):
    fop=open(firstorderpath,'r')
    fox=[]
    foy=[]
    folen=[]
    fowid=[]
    foid=[]
    for line in fop:
        #if line[0]!='#':
        linesplit = line.split() 
        fox.append(float( linesplit[1][0:-1] ))  ### [0:-1] strips off the comma. 
        foy.append(float(linesplit[2][0:-1]))
        folen.append(float(linesplit[3][0:-1])) 
        fowid.append(float(linesplit[-1].split('{')[-1].split('}')[0]))  ### python is weird.
            #fowid.append(float(line[26]))
        #elif line[0]=='#':
        #    ll=len(line)
         #   indend=ll-2
         #   ind1=indend-6
         #   for i in range(7):
         #       if line[ind1+i]!=' ' and line[ind1+i]!='{':
         #           indstart=ind1+i
         #           break
            #print line[indstart:85]
         #   if indstart<indend:
         #       foid.append(int(line[indstart:indend]))
         #   else:
           #     foid.append(int(line[indend]))
    #fop.close()
    return (fox,foy,folen,fowid,foid)


def go (linelistfile=" ",save_temp=True,recover_temp=False,show_dispersed=False):
    if linelistfile==" ":
        allDirectoryFiles=os.listdir('.')
        for files in allDirectoryFiles:
            if files[0:3]=='Par':
                llpts=files.split('_')
                linelistfile='linelist/'+llpts[0]+'lines.dat'
                break
    if linelistfile==" " or os.path.exists(linelistfile)==0:
        print "Invalid path to line list file: %s" % (linelistfile)
        return 0
    else:
        print "Found line list file %s" % (linelistfile)

    ### CREAT TRACE REGION FILES ###
    trace102 = open('G102_trace.reg', 'w')
    trace102.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace102.write('wcs;\n')
    trace102.write('box(9950,0,3100,1,1.62844e-12)\n')
    trace102.close()
    trace141 = open('G141_trace.reg', 'w')
    trace141.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace141.write('wcs;\n')
    trace141.write('box(14142.49,0,5500,1,0)\n')
    trace141.close()

#    # check for trace region files
#    if os.path.isfile('G102_trace.reg')==0 or os.path.isfile('G141_trace.reg')==0:
#        print '\nMissing one or both trace region files.'
#        print 'Expecting them in Spectra/ directory.'
#        print 'Exiting'
#        return 0

    lam_Halpha=6563.0
    lam_Hbeta=4861.0
    lam_Oiii_1=4959.0
    lam_Oiii_2=5007.0
    lam_Oii=3727.0
    lam_Sii=6724.0
    lam_Siii_1=9069.0
    lam_Siii_2=9532.0
    lam_Lya=1216.0
    lam_He=10830.0
    lam_Fe=12600.0
    lam_Pag=10940.0
    lam_Pab=12810.0
    
    suplines=[lam_Lya,lam_Oii,lam_Hbeta,lam_Oiii_1,lam_Oiii_2,lam_Halpha,lam_Sii,lam_Siii_1,lam_Siii_2,lam_He,lam_Pag,lam_Fe,lam_Pab]
    
    parnos,grism,objid,wavelen,npix,ston,flag,nstart,setzs,flagcont,comment=readll(linelistfile,recover_temp)
    
    parts=linelistfile.split('.')
    linelistoutfile=parts[0]+'_with_redshifts.'+parts[1]
    linelistoutfile_cont=parts[0]+'_with_redshifts_contam.'+parts[1] #MR
    # MB: avoid overwriting temp file by accident
    lltemp=linelistoutfile[0:-3]+'incomplete.dat'
    lctemp=linelistoutfile_cont[0:-3]+'incomplete.dat'
    if os.path.isfile(lltemp):
        print '\nTemp files: \n  %s,\n  %s \nalready exist and you are not set to recover them.\n' % (lltemp, lctemp)
        ask = raw_input('Overwrite? [y/n] ')
        if ask.lower() == 'n':
            print 'Exiting.'
            return 0
        

    g102zeroarr=[]
    g102firstarr=[]
    g102zeroordreg="../DATA/DIRECT_GRISM/G102_0th.reg"
    g102firstordreg="../DATA/DIRECT_GRISM/G102_1st.reg"
    if os.path.exists(g102zeroordreg)==1:
        g102zeroarr=getzeroorders(g102zeroordreg,g='G102')
        g102firstarr=getfirstorders(g102firstordreg)

    g141zeroarr=[]
    g141firstarr=[]
    g141zeroordreg="../DATA/DIRECT_GRISM/G141_0th.reg"
    g141firstordreg="../DATA/DIRECT_GRISM/G141_1st.reg"
    if os.path.exists(g141zeroordreg)==1:
        g141zeroarr=getzeroorders(g141zeroordreg)
        g141firstarr=getfirstorders(g141firstordreg)
        
    if len(g102zeroarr)>0:
        g102zerox=g102zeroarr[0]
        g102zeroy=g102zeroarr[1]
        g102zeroid=g102zeroarr[2]
        g102firstx=g102firstarr[0]
        g102firsty=g102firstarr[1]
        g102firstlen=g102firstarr[2]
        g102firstwid=g102firstarr[3]
        g102firstid=g102firstarr[4]
    else:
        g102zerox=[]
        g102zeroy=[]
        g102zeroid=[]
        g102firstx=[]
        g102firsty=[]
        g102firstlen=[]
        g102firstwid=[]
        g102firstid=[]
        
    if len(g141zeroarr)>0:
        g141zerox=g141zeroarr[0]
        g141zeroy=g141zeroarr[1]
        g141zeroid=g141zeroarr[2]
        g141firstx=g141firstarr[0]
        g141firsty=g141firstarr[1]
        g141firstlen=g141firstarr[2]
        g141firstwid=g141firstarr[3]
        g141firstid=g141firstarr[4]
    else:
        g141zerox=[]
        g141zeroy=[]
        g141zeroid=[]
        g141firstx=[]
        g141firsty=[]
        g141firstlen=[]
        g141firstwid=[]
        g141firstid=[]
        

    zset=0
    option=' '
    #setzs=[] # moved to readll subroutine by NR for version 1.0
    #flagcont=[] #MR # moved to readll subroutine by NR for version 1.0
    #comment=[] #MR # moved to readll subroutine by NR for version 1.0

    #for i in range(len(parnos)): # moved to readll subroutine by NR for version 1.0
    #    setzs.append(0.0) # moved to readll subroutine by NR for version 1.0
    #    flagcont.append(0) #MR # moved to readll subroutine by NR for version 1.0
    #    comment.append('') #MR # moved to readll subroutine by NR for version 1.0

    i=nstart
    #print i
    #for i in range(len(parnos)):
    if recover_temp:
        showDirectNEW(1,0)
        if show_dispersed:  # MB
            showDispersed(1,0)
    progress=0.0
    while i<len(parnos):
        specname='Par' + parnos[i] + '_BEAM_' + str(objid[i]) + 'A.dat'
        specnameg102='Par' + parnos[i] + '_G102_BEAM_' + str(objid[i]) + 'A.dat'
        specnameg141='Par' + parnos[i] + '_G141_BEAM_' + str(objid[i]) + 'A.dat'
        plotTitle='Par' + parnos[i] + '\_BEAM\_' + str(objid[i])
        if os.path.exists(specnameg102)==0:
            specname=specnameg141
        cname102='linelist/Par' + parnos[i] + '_G102_BEAM_' + str(objid[i]) + 'A.contfit.dat'
        cname141='linelist/Par' + parnos[i] + '_G141_BEAM_' + str(objid[i]) + 'A.contfit.dat'
        
        lamline=wavelen[i]
        progress=float(i)/float(len(parnos))*100.0
        print "Progress: %.1f percent" % (progress)
        if i<0:
            i=0

        if i==nstart or objid[i]!=objid[i-1]:
            zset=0
            speclam=[]
            specval=[]
            specunc=[]
            speccon=[]
            speczer=[]
            if specname==specnameg141:
                speclam,specval,specunc,speccon,speczer=readspec(specname,isG141=True)
            else:
                speclam,specval,specunc,speccon,speczer=readspec(specname,isG141=False)
            continlam=[]
            continflam=[]

            if os.path.exists(cname102)==1:
                c102in=open(cname102,'r')
                for line in c102in:
                    entries=line.split()
                    continlam.append(float(entries[0]))
                    continflam.append(float(entries[1]))
                c102in.close()
            if os.path.exists(cname141)==1:
                c141in=open(cname141,'r')
                for line in c141in:
                    entries=line.split()
                    continlam.append(float(entries[0]))
                    continflam.append(float(entries[1]))
                c141in.close()
            spec_lam=np.array(speclam)
            spec_val=np.array(specval)
            spec_unc=np.array(specunc)
            spec_con=np.array(speccon)
            spec_zer=np.array(speczer)
            zero_ord = spec_val * spec_zer
            if len(continlam)>0:
                cont_val=scipy.interpolate.spline(np.array(continlam),np.array(continflam),spec_lam)
            
            if len(g102zerox)>0:
                show2dNEW('G102',parnos[i],int(objid[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear')
                
            if len(g141zerox)>0:
                show2dNEW('G141',parnos[i],int(objid[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear')
                showDirectNEW(objid[i],i)
                if show_dispersed:  # MB
                    showDispersed(objid[i],i)
        
        next=0

        xmin=spec_lam.min()-200.0
        xmax=spec_lam.max()+200.0
        ymin=spec_val.min()
        ymax=1.5*spec_val.max()
    
        if i==nstart or zset==0:
            zguess=(lamline/lam_Halpha)-1
            slineslam=[]
            slinesvals=[]
        
        linex=np.array([lamline,lamline])
        liney=np.array([ymin,ymax])
        while (next==0):
            linevecs=ret_z_lines(suplines,zguess,xmin,xmax,ymin,ymax)
            print "Guessing z = %f" % (zguess)
            supx=np.array(linevecs[0])
            supy=np.array(linevecs[1])
            
            linptx=np.array([lamline])
            linpty=np.array([(ymin+ymax)/2.0])
            
            plt.ion()
            plt.figure(1,figsize=(11,8))
            plt.clf()
            plt.subplot(211)
            if len(continlam)>0:
                cplt=plt.plot(np.array(continlam),np.array(continflam),'b--')
            specplt=plt.plot(spec_lam, spec_val, 'k',spec_lam,spec_con,'r',spec_lam,zero_ord,'m',supx,supy,'g:')
            plt.setp(specplt, ls = 'steps')
            linepts=plt.scatter(linptx,linpty,s=80,c='g',marker='o', cmap=None, norm=None,vmin=None, vmax=None, alpha=1, linewidths=None,verts=None)
            plt.axvline(x=11500.,c='c',linestyle=':')
            plt.xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
            plt.ylabel(r'F$_\lambda$ ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
            plt.xlim([xmin, xmax])
            plt.ylim([ymin, ymax])
            plt.title(plotTitle)
            if len(slineslam)>0:
                linesplt=plt.plot(slines_lam,slines_val,'g:')
            plt.subplot(212)
            s2n=(spec_val-cont_val)/spec_unc
            s2n_lam=spec_lam
            mask=np.logical_and(s2n>-10000., s2n<10000.)
            s2n=s2n[mask]
            s2n_lam=s2n_lam[mask]
            plt.plot(s2n_lam,s2n,'k-',linestyle='steps')
            plt.axhline(y=1.73,c='r')
            plt.xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
            plt.ylabel(r'S/N',size='xx-large')
            plt.xlim([xmin, xmax])
            
            plt.draw()
            plt.draw()
            if (option!='s' or objid[i]!=objid[i-1]):
                print "Enter option: s=skip spectrum b=back z=set redshift zc=set redshift with contamination r=reject line rc=reject line and mark contaminated c=leave comment ha=Halpha hb=Hbeta o2=OII o31=OIII4959 o32=OIII5007 s2=SII s31=SIII9069 s32=SIII9532 la=Lyman alpha mb=shift line blueward mr=shift line redward q=quit lin=linear z-scale log=logarithmic zscale zs102=z1,z2 comma-separated range for G102 zscale zs141=z1,z2 comma-separated range for G141 zscale dc=recenter direct images dr=reload direct image reg files" #MR
                option = raw_input(">")
            if option=='r':
                next=1
            elif option=='z' and zset==0:
                next=1
                zset=1
                flagcont[i]=1 #MR
                setzs[i]=zguess
                zprev=zguess
            elif option=='z' and zset==1:
                flagcont[i]=1 #MR
                setzs[i]=zprev
                next=1
            #MR Added 10-25-11
            #elif option=='c' and zset==0:
	    elif option=='c': # Changed by NR version 0.2.5
                print "Enter your comment here:"
                comment[i]=raw_input(">")
            elif option=='zc' and zset==0:
                next=1
                zset=1
                flagcont[i]=2
                setzs[i]=zguess
                zprev=zguess
            elif option=='zc' and zset==1:
                flagcont[i]=2
                setzs[i]=zprev
                next=1
            elif option=='rc':
                next=1
                flagcont[i]=3
            elif option=='mbl':
                print "Shifting line 100Angs blueward"
                #wavelen[i]=wavelen[i]-10.0
                lamline=lamline-100.0
                wavelen[i]=lamline
            elif option=='mrl':
                print "Shifting line 100Angs redward"
                #wavelen[i]=wavelen[i]+10.0
                lamline=lamline+100.0
                wavelen[i]=lamline
            # End MR added
            elif option=='ha':
                zguess=(lamline/lam_Halpha)-1
            elif option=='hb':
                zguess=(lamline/lam_Hbeta)-1
            elif option=='o2':
                zguess=(lamline/lam_Oii)-1
            elif option=='o31':
                zguess=(lamline/lam_Oiii_1)-1
            elif option=='o32':
                zguess=(lamline/lam_Oiii_2)-1
            elif option=='s2':
                zguess=(lamline/lam_Sii)-1
            elif option=='s31':
                zguess=(lamline/lam_Siii_1)-1
            elif option=='s32':
                zguess=(lamline/lam_Siii_2)-1
            elif option=='la':
                zguess=(lamline/lam_Lya)-1
            elif option=='b':
                print "Going back."
                i=i-2
                next=1
            elif option=='s':
                print "Skipping ahead to next spectrum."
                next=1
            elif option=='mb':
                print "Shifting line 10Angs blueward"
                #wavelen[i]=wavelen[i]-10.0
                lamline=lamline-10.0
                wavelen[i]=lamline
            elif option=='mr':
                print "Shifting line 10Angs redward"
                #wavelen[i]=wavelen[i]+10.0
                lamline=lamline+10.0
                wavelen[i]=lamline
            elif option=='q':
                print "Quitting without saving"
                return 0
            elif option=='lin':
                if len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear')
                if len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear')
            elif option=='log':
                if len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'log')
                if len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'log')
            elif len(option)>6 and option[0:6]=='zs102=':
                vals=option[6:]
                zran=vals.split(',')
                if len(zran)!=2 or isFloat(zran[0])==False or isFloat(zran[1])==False:
                    print "Invalid zrange."
                elif len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear',zran1=float(zran[0]),zran2=float(zran[1]))
            elif len(option)>6 and option[0:6]=='zs141=':
                vals=option[6:]
                zran=vals.split(',')
                if len(zran)!=2 or isFloat(zran[0])==False or isFloat(zran[1])==False:
                    print "Invalid zrange."
                elif len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear',zran1=float(zran[0]),zran2=float(zran[1]))
            elif option=='dc':
                showDirectNEW(objid[i],i)
                if show_dispersed:  # MB
                    showDispersed(objid[i],i)
            elif option=='dr':
                reloadReg()
            else:
                print "Invalid entry.  Try again."
            print "OK"
        if save_temp:
            lltemp=linelistoutfile[0:-3]+'incomplete.dat'
            lctemp=linelistoutfile_cont[0:-3]+'incomplete.dat'
            printLLouttemp(lltemp,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
            printLCouttemp(lctemp,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
            
        i=i+1
    printLLout(linelistoutfile,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
    printLCout(linelistoutfile_cont,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
    # Clean up temp files
    if save_temp:
        if os.path.exists(lltemp)==1:
            os.unlink(lltemp)
        if os.path.exists(lctemp)==1:
            os.unlink(lctemp)
    if os.path.exists('./tempcoo.dat')==1:
        os.unlink('./tempcoo.dat')
    if os.path.exists('./temp_zero_coords.coo')==1:
        os.unlink('./temp_zero_coords.coo')
    if os.path.exists('./temp110.fits')==1:
        os.unlink('./temp110.fits')
    if os.path.exists('./temp140.fits')==1:
        os.unlink('./temp140.fits')
    if os.path.exists('./temp160.fits')==1:
        os.unlink('./temp160.fits')
    if os.path.exists('./temp_zero_coords.reg')==1:
        os.unlink('./temp_zero_coords.reg')
    

def printLLout(linelistoutfile,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment): #Edited for version 1.0 by NR
    if os.path.exists(linelistoutfile)==1:
        os.unlink(linelistoutfile)
    llout=open(linelistoutfile,'w')
    for i in range(len(parnos)):
        if setzs[i]>0.0:
            print >> llout, "%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%f\t%s" % (parnos[i],grism[i],objid[i],wavelen[i],npix[i],ston[i],flag[i],flagcont[i],setzs[i],comment[i]) #MR
    llout.close()
  

def printLCout(linelistoutfile_cont,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment):
    if os.path.exists(linelistoutfile_cont)==1:
        os.unlink(linelistoutfile_cont)

    llout=open(linelistoutfile_cont,'w')
    for i in range(len(parnos)):
        if flagcont[i]>2:
            print >> llout, "%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%f\t%s" % (parnos[i],grism[i],objid[i],wavelen[i],npix[i],ston[i],flag[i],flagcont[i],setzs[i],comment[i])
    llout.close()
    # End MR added

def printLLouttemp(linelistoutfile,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment): #Edited for version 1.0 by NR
    if os.path.exists(linelistoutfile)==1:
        os.unlink(linelistoutfile)
    llout=open(linelistoutfile,'w')
    for i in range(len(parnos)):
        if setzs[i]>0.0:
            print >> llout, "%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%f\t%i\t%s" % (parnos[i],grism[i],objid[i],wavelen[i],npix[i],ston[i],flag[i],flagcont[i],setzs[i],i,comment[i]) #MR
    llout.close()
  

def printLCouttemp(linelistoutfile_cont,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment):
    if os.path.exists(linelistoutfile_cont)==1:
        os.unlink(linelistoutfile_cont)

    llout=open(linelistoutfile_cont,'w')
    for i in range(len(parnos)):
        if flagcont[i]>2:
            print >> llout, "%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%f\t%i\t%s" % (parnos[i],grism[i],objid[i],wavelen[i],npix[i],ston[i],flag[i],flagcont[i],setzs[i],i,comment[i])
    llout.close()
    # End MR added

def show2dNEW (grism,parno,obid,firstx,firsty,firstlen,firstwid,firstid,zerox,zeroy,zeroid,trans,zran1=-0.2,zran2=0.75):
# In version 1.0, will first look for wavelength-calibrated stamps in the G1??_DRIZZLE directories; failing this, will default to old stamps
    dims=()
    zrad=10.0
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'
    path2dl=par_root_dir + grism + '_DRIZZLE/aXeWFC3_' +grism + '_mef_ID'+str(obid)+'.fits'
    if os.path.exists(path2dl)==1:
        path2d=path2dl
    else:
        path2d=par_root_dir+'Stamps/Par'+ str(parno)+'_'+grism+'_BEAM_'+str(obid)+'A.fits'

    if grism=='G102':
        frameno='1'
    elif grism=='G141':
        frameno='2'
    if os.path.exists(path2d)==1:
        infits=pyfits.open(path2d)
        ### changing to read in 1st data extension ###
        #darr=infits[-1].data
        hdr=infits[1].header
        darr=infits[1].data
        dims=darr.shape
        infits.close()
    elif os.path.exists(path2d)==0:
        print "%s stamp not found." % (grism)
        return False

    ### USING THE DRIZZLE TRANSFORMATIONS TO GET ZEROTH ORDERS ###
    _cx = np.array([xcoo for xcoo in zerox]) - hdr['BB0X'] - 1
    _cy = np.array([ycoo for ycoo in zeroy]) - hdr['BB0Y'] - 1
    cx = hdr['D001OUXC'] + (hdr['DRZ00'] + hdr['DRZ01']*(_cx-hdr['D001INXC']) + hdr['DRZ02']*(_cy-hdr['D001INYC']))
    cy = hdr['D001OUYC'] + (hdr['DRZ10'] + hdr['DRZ11']*(_cx-hdr['D001INXC']) + hdr['DRZ12']*(_cy-hdr['D001INYC']))
    # convert to (Angs,arcsec) coords
    cx = (cx - hdr['CRPIX1'])*hdr['CDELT1'] + hdr['CRVAL1']
    cy = (cy - hdr['CRPIX2'])*hdr['CDELT2'] + hdr['CRVAL2']
    rad = 5 * hdr['CDELT1']
    outcoo=par_root_dir+"Spectra/temp_zero_coords.reg"
    if os.path.exists(outcoo)==1:
        os.unlink(outcoo)
    f = open(outcoo, 'w')
    f.write('wcs;\n')
    for j in range(len(zerox)):
        f.write('circle(%.2f,%.4f,%.1f) # color=red text={%s}\n' % (cx[j],cy[j],rad,zeroid[j]))
    f.close()

    ### THIS WAS ALL FOR THE OLD TRANSFORMATION ###
#    matchind=0
#    i=0
#    for fid in firstid:
#        if fid==obid:
#            matchind=i
#            break
#        i=i+1
#    xmin=firstx[matchind]-firstlen[matchind]/2.0
#    xmax=firstx[matchind]+firstlen[matchind]/2.0
#    ymin=firsty[matchind]-firstwid[matchind]/2.0
#    ymax=firsty[matchind]+firstwid[matchind]/2.0
#    
#    numzer=0
#    if len(dims)>0:
#        xdim=float(max(dims))
#        ydim=float(min(dims))
#    outcoo=par_root_dir+"Spectra/temp_zero_coords.reg"
#    if os.path.exists(outcoo)==1:
#        os.unlink(outcoo)
#    coordout=open(outcoo,'w')
#    print >>coordout, "image"
#    for j in range(len(zerox)):
#        # MB: use dims of stamp, rather than size of 1st order in region file
#        if zerox[j] >= (firstx[matchind] - xdim/2.) and \
#           zerox[j] <= (firstx[matchind] + xdim/2.) and \
#           zeroy[j] >= (firsty[matchind] - ydim/2.) and \
#           zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G102':
##    	if zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G102':
#    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+212./2.0-13,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+3.6,zeroid[j])
#
#        elif zerox[j] >= (firstx[matchind] - xdim/2.) and \
#             zerox[j] <= (firstx[matchind] + xdim/2.) and \
#             zeroy[j] >= (firsty[matchind] - ydim/2.) and \
#             zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G141':
##    	elif  zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G141':
#    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+184./2.0,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+0.6,zeroid[j])
#    numzer=numzer+1
#    coordout.close()
    
    if trans=='log':
        zscale='log'
    else:
        zscale='linear'
    cmd='xpaset -p ds9 frame '+frameno
    os.system(cmd)
    cmd='xpaset -p ds9 file '+path2d
    os.system(cmd)
    cmd='xpaset -p ds9 scale limits '+str(zran1)+' '+str(zran2)
    os.system(cmd)
    cmd='xpaset -p ds9 scale '+zscale
    os.system(cmd)
    #cmd='xpaset -p ds9 zoom to fit'
    #os.system(cmd)
    cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/temp_zero_coords.reg'
    os.system(cmd)
    # MR
    if frameno=='1':
        cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G102_trace.reg'
    if frameno=='2':
        cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G141_trace.reg'
    os.system(cmd)
    
def showDirectNEW(obid,lineno=-1):
    obid=int(obid)
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'

    path2direct=par_root_dir+'DATA/DIRECT_GRISM/'
    #path110=path2direct+'F110W_rot_drz.fits'
    #path140=path2direct+'F140W_rot_drz.fits'
    #path160=path2direct+'F160W_rot_drz.fits'
    path110=path2direct+'F110W_drz.fits'
    path140=path2direct+'F140W_drz.fits'
    path160=path2direct+'F160W_drz.fits'
    path140cat=path2direct+'fin_F140.cat'
    path160cat=path2direct+'fin_F160.cat'
    if os.path.exists(path110)==0 and os.path.exists(path140)==0 and os.path.exists(path160)==0:
        print "No Direct Images Found."
        return 0
    if os.path.exists(path140cat)==1:
        infHcat=open(path140cat,'r')
    elif os.path.exists(path160cat)==1:
        infHcat=open(path160cat,'r')
    else:
        return 0
    xcen,ycen=-1,-1
    for line in infHcat:
        if line[0]!='#':
            entries=line.split()
            if int(entries[1])==obid:
                xcenter,ycenter=float(entries[7]),float(entries[8])
                hexcoo=[entries[7],entries[8]]
    infHcat.close()
    if lineno!=0 and os.path.exists(path110)==1:
        panDirect(hexcoo[0],hexcoo[1],grism='G102')
    if lineno!=0:
        panDirect(hexcoo[0],hexcoo[1])
        return 0
    if os.path.exists(path110)==1:
        cmd='xpaset -p ds9 frame 3'
        os.system(cmd)
        cmd='xpaset -p ds9 file '+path110
        os.system(cmd)
        ### using F110_drz.reg with F110W_drz.fits
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg'
        os.system(cmd)
        cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
        os.system(cmd)
    if os.path.exists(path140)==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 file '+path140
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+ 'DATA/DIRECT_GRISM/F140_drz.reg'
        os.system(cmd)
        cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
        os.system(cmd)
    elif os.path.exists(path160)==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 file '+path160
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg'
        os.system(cmd)
        cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
        os.system(cmd)
    

def showDispersed(obid,lineno=-1):  # MB
    obid=int(obid)
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'

    path2dispersed=par_root_dir+'DATA/DIRECT_GRISM/'
    ### Using G102.fits instead of G102_drz.fits ###
    path102=path2dispersed+'G102.fits'
    path141=path2dispersed+'G141.fits'
    path102_0reg = os.path.join(path2dispersed, 'G102_0th.reg')
    path102_1reg = os.path.join(path2dispersed, 'G102_1st.reg')
    path141_0reg = os.path.join(path2dispersed, 'G141_0th.reg')
    path141_1reg = os.path.join(path2dispersed, 'G141_1st.reg')
    if os.path.exists(path102)==0 and os.path.exists(path141)==0:
        print "No Grism Images Found."
        return 0
    # get center of 1st order
    ### Using same syntax as getzeroorders for consistency ###
    if os.path.exists(path102_1reg)==1:
        reg102=open(path102_1reg,'r')
        x102,y102=-1,-1
        for line in reg102:
            # using same syntax as getzeroorders
            linesplit = line.split()
            textid = int(re.search('\d+',linesplit[-1]).group(0))
            if textid==obid:
                x102 = float(linesplit[1].split(',')[0])
                y102 = float(linesplit[2].split(',')[0])
        reg102.close()
    else:
        return 0
    if os.path.exists(path141_1reg)==1:
        reg141=open(path141_1reg,'r')
        x141,y141=-1,-1
        for line in reg141:
            linesplit = line.split()
            textid = int(re.search('\d+',linesplit[-1]).group(0))
            if textid==obid:
                x141 = float(linesplit[1].split(',')[0])
                y141 = float(linesplit[2].split(',')[0])
        reg141.close()
    else:
        return 0
    
    if lineno!=0 and os.path.exists(path102)==1:
        panDispersed(x102,y102,grism='G102')
    if lineno!=0:
        panDispersed(x141,y141)
        return 0
    if os.path.exists(path102)==1:
        cmd='xpaset -p ds9 frame 5'
        os.system(cmd)
        cmd='xpaset -p ds9 file '+path102
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+path102_0reg
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+path102_1reg
        os.system(cmd)
        cmd='xpaset -p ds9 pan to %f %f image' % (x102,y102)
        os.system(cmd)
    if os.path.exists(path141)==1:
        cmd='xpaset -p ds9 frame 6'
        os.system(cmd)
        cmd='xpaset -p ds9 file '+path141
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+path141_0reg
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+path141_1reg
        os.system(cmd)
        cmd='xpaset -p ds9 pan to %f %f image' % (x141, y141)
        os.system(cmd)


def createAltGrismRegion(grism):
    workingdir=os.getcwd()
    par = os.path.dirname(workingdir)
    cat110 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F110.cat')
    cat140 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F140.cat')
    cat160 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F160.cat')

    if grism == 'G102':
        if os.path.exists(cat110) == 1:
            cat = np.genfromtxt(cat110, dtype=[('num',int),('a_img',float),
                               ('mag',float)], usecols=(1,4,12))
        else:
            print cat110
            return 0
    if grism == 'G141':
        if os.path.exists(cat140) == 1:
            cat = np.genfromtxt(cat140, dtype=[('num',int),('a_img',float),
                                 ('mag',float)], usecols=(1,4,12))
        elif os.path.exists(cat160) == 1:
            cat = np.genfromtxt(cat160, dtype=[('num',int),('a_img',float),
                                 ('mag',float)], usecols=(1,4,12))
        else:
            print 'nope2'
            return 0

    f = open(os.path.join(workingdir,grism+'_temp.reg'), 'w')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')
    cenx = np.zeros(cat['num'].shape[0], dtype=float)
    ceny = np.zeros(cat['num'].shape[0], dtype=float)
    for i in range(cat['num'].shape[0]):
        objid = cat['num'][i]
        gfile = os.path.join(par,'%s_DRIZZLE/aXeWFC3_%s_mef_ID%i.fits'%(grism,grism,objid))
        if os.path.isfile(gfile):
            hdr = pyfits.getheader(gfile, 1)
            # coords of bounding box
            boxx = np.array([hdr['bb0x'], hdr['bb1x']])
            boxy = np.array([hdr['bb0y'], hdr['bb1y']])
            cenx[i] = boxx[0] + (boxx[1] - boxx[0]) / 2.
            ceny[i] = boxy[0] + (boxy[1] - boxy[0]) / 2.
            slitwidth = hdr['slitwidt']
            sx = 184
            sy = 8
#            sy = slitwidth * cat['a_img'][i]
            f.write('box(%f,%f,%i,%i,0) # text={%i}\n' % (cenx[i],ceny[i],sx,sy,objid))
    f.close()
    return cenx,ceny


def panDirect(ra,dec,grism='G141'):
    # Pan to coords in frame
    if grism=='G141':
        fno='4'
    else:
        fno='3'
    cmd='xpaset -p ds9 frame ' + fno
    os.system(cmd)
    cmd='xpaset -p ds9 pan to '+ra+' '+dec+' fk5'
    os.system(cmd)


def panDispersed(xx,yy,grism='G141'):  # MB
    # Pan to coords in frame
    if grism=='G141':
        fno='6'
    else:
        fno='5'
    cmd='xpaset -p ds9 frame ' + fno
    os.system(cmd)
    cmd='xpaset -p ds9 pan to %f %f image' % (xx,yy)
    os.system(cmd)


def reloadReg():
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'

    if os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F110.reg')==1:
        cmd='xpaset -p ds9 frame 3'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110.reg'
        os.system(cmd)
    if os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F160.reg')==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160.reg'
        os.system(cmd)
    elif os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F140.reg')==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F140.reg'
        os.system(cmd)

def readspec(bpath,isG141=False):
    darray=np.genfromtxt(bpath,dtype=np.float)
    l=darray[0:,0]
    fl=darray[0:,1]
    efl=darray[0:,2]
    cfl=darray[0:,3]
    zfl=darray[0:,4]
    maskNaN_fl=np.isnan(fl)
    maskNaN_efl=np.isnan(efl)
    maskNaN=np.logical_or(maskNaN_fl,maskNaN_efl)
    notNaN=np.logical_not(maskNaN)
    if isG141:
        masklam=np.logical_and(l>11000.,l<17000.)
        notNaN=np.logical_and(masklam,notNaN)
    l,fl,efl,cfl,zfl=l[notNaN],fl[notNaN],efl[notNaN],cfl[notNaN],zfl[notNaN]
    return l,fl,efl,cfl,zfl

def readll(linelistfile,recover_temp):
    parnos=[]
    grism=[]
    objid=[]
    wavelen=[]
    npix=[]
    ston=[]
    flag=[]
    
    llin=open(linelistfile,'r')
    for line in llin:
        entries=line.split()
        parnos.append(entries[0])
        grism.append(entries[1])
        objid.append(int(entries[2]))
        wavelen.append(float(entries[3]))
        npix.append(entries[4])
        ston.append(entries[5])
        flag.append(entries[6])
    llin.close()
    setzs=[]
    flagcont=[] #MR
    comment=[] #MR

    for i in range(len(parnos)):
        setzs.append(0.0)
        flagcont.append(0) #MR
        comment.append('') #MR

    if recover_temp==False:
        return np.array(parnos),np.array(grism),np.array(objid),np.array(wavelen),np.array(npix),np.array(ston),np.array(flag),0,np.array(setzs),np.array(flagcont),np.array(comment, dtype='S100')
    else:
        nstart1,nstart2=0,0
        parts=linelistfile.split('.')
        linelistoutfile=parts[0]+'_with_redshifts.dat'
        lltemp=linelistoutfile[0:-3]+'incomplete.dat'
        lctemp=linelistoutfile[0:-4]+'_contam.incomplete.dat'
        nstart=0
        print "Attempting to recover progress from files %s and %s" % (lltemp,lctemp)
        if os.path.exists(lltemp)==0 or os.path.exists(lctemp)==0:
            return np.array(parnos),np.array(grism),np.array(objid),np.array(wavelen),np.array(npix),np.array(ston),np.array(flag),nstart,np.array(setzs),np.array(flagcont),np.array(comment, dtype='S100')
        # INCOMPLETE
        # Need to figure out how to merge in saved temp files!
        inftemp=open(lltemp,'r')
        for line in inftemp:
            if line[0]!='#':
                entries=line[0:-1].split('\t') # should delete the \n character from the end.
                line_ind=int(entries[9])
                nstart1=line_ind
                wavelen[line_ind]=float(entries[3])
                flagcont[line_ind]=int(entries[7])
                setzs[line_ind]=float(entries[8])
                if len(entries)==11 and len(entries[10])>1:
                    comment[line_ind]=entries[10]
        inftemp.close()
        inftemp=open(lctemp,'r')
        for line in inftemp:
            if line[0]!='#':
                entries=line[0:-1].split('\t')# should delete the \n character from the end.
                line_ind=int(entries[9])
                nstart2=line_ind
                wavelen[line_ind]=float(entries[3])
                flagcont[line_ind]=int(entries[7])
                setzs[line_ind]=float(entries[8])
                if len(entries)==11 and len(entries[10])>1:
                    comment[line_ind]=entries[10] 
        inftemp.close()
        nstart=max([nstart1,nstart2])+1 # starting from the next line candidate
        #print nstart
        return np.array(parnos),np.array(grism),np.array(objid),np.array(wavelen),np.array(npix),np.array(ston),np.array(flag),nstart,np.array(setzs),np.array(flagcont),np.array(comment, dtype='S100')
