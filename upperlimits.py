import numpy as np
import scipy
from scipy.integrate import quad
import os
import distutils
import math

def linlsq(xs,ys,sigs):
    sumsig=0.0
    sumx=0.0
    sumy=0.0
    sumx2=0.0
    sumxy=0.0
    for j in range(len(sigs)):
        sumsig=sumsig+(1.0/sigs[j]**2)
        sumx=sumx+(xs[j]/sigs[j]**2)
        sumy=sumy+(ys[j]/sigs[j]**2)
        sumx2=sumx2+(xs[j]**2/sigs[j]**2)
        sumxy=sumxy+(xs[j]*ys[j]/sigs[j]**2)
    fmat=np.array([[sumsig,sumx],[sumx,sumx2]])
    fans=np.array([sumy,sumxy])
    if fmat[0][0]*fmat[1][1]-fmat[1][0]*fmat[0][1]==0.0:
        params=[0.0,0.0]
    else:
        #print fmat[0][0]*fmat[1][1]-fmat[1][0]*fmat[0][1]
        finv=np.linalg.inv(fmat)
        fanst=np.transpose(fans)
        params=np.dot(finv,fanst) # Line fit with form params[0]+params[1]*x
    return params

def getUL(l0line,l0ref,zref,FWHMref,specpath):
    npix=15
    sigref=FWHMref/2.35
    infUL=open(specpath,'r')
    speclam=[]
    specval=[]
    specunc=[]
    for lines in infUL:
        ents=lines.split()
        if ents[1]!='nan' and ents[1]!='NaN' and ents[1]!='NAN' and ents[1]!='-nan' and ents[1]!='-NaN' and ents[1]!='-NAN' and ents[2]!='nan' and ents[2]!='NaN' and ents[2]!='NAN' and ents[2]!='-nan' and ents[2]!='-NaN' and ents[2]!='-NAN':
            speclam.append(float(ents[0]))
            specval.append(float(ents[1]))
            specunc.append(float(ents[2]))
    infUL.close()
    if l0line*(1.0+zref)<min(speclam) or l0line*(1.0+zref)>max(speclam):
        return 0
    if l0line<l0ref:
        lamUL=(l0line-2.5*sigref)*(1.0+zref)
        lims=[min(speclam),lamUL]
    elif l0line>l0ref:
        lamUL=(l0line+2.5*sigref)*(1.0+zref)
        lims=[lamUL,max(speclam)]
    speclamclipped=[]
    specvalclipped=[]
    specuncclipped=[]
    for i in range(len(speclam)):
        if speclam[i]>lims[0] and speclam[i]<lims[1]:
            speclamclipped.append(speclam[i])
            specvalclipped.append(specval[i])
            specuncclipped.append(specunc[i])

    nclipped=len(speclamclipped)
    if nclipped>npix and l0line<l0ref:
        speclamclipped=speclamclipped[nclipped-npix-1:]
        specvalclipped=specvalclipped[nclipped-npix-1:]
        specuncclipped=specuncclipped[nclipped-npix-1:]
    elif nclipped>npix and l0line>l0ref:
        speclamclipped=speclamclipped[0:npix]
        specvalclipped=specvalclipped[0:npix]
        specuncclipped=specuncclipped[0:npix]
    if len(speclamclipped)<7:
        return 0.0
    fitparams=linlsq(speclamclipped,specvalclipped,specuncclipped)
    #fitparams=np.polyfit(speclamclipped,specvalclipped)
    if fitparams[0]==0.0 and fitparams[1]==0.0:
        return 0.0
    pixvar=0.0
    for i in range(len(speclamclipped)):
        pixvar=pixvar+(specvalclipped[i]-(fitparams[0]+fitparams[1]*speclamclipped[i]))**2
    deltalam=speclamclipped[1]-speclamclipped[0]
    #fsig=(2.35482*sigref/deltalam)**0.5*(pixvar/(float(len(speclamclipped))-2.0))**0.5*1.564206
    fsig=(2.35482*sigref*(1.0+zref)*deltalam)**0.5*(pixvar/(float(len(speclamclipped))-2.0))**0.5*1.564206
    return fsig

def gauss(x):
    return (2.*math.pi)**(-0.5)*math.exp(-x**2/2.0)

def getULerrspec(l0line,zref,FWHMref,specpath):
    #npix=15
    #sigref=FWHMref/2.35
    infUL=open(specpath,'r')
    speclam=[]
    specval=[]
    specunc=[]
    for lines in infUL:
        ents=lines.split()
        if ents[1]!='nan' and ents[1]!='NaN' and ents[1]!='NAN' and ents[1]!='-nan' and ents[1]!='-NaN' and ents[1]!='-NAN' and ents[2]!='nan' and ents[2]!='NaN' and ents[2]!='NAN' and ents[2]!='-nan' and ents[2]!='-NaN' and ents[2]!='-NAN' and float(ents[0])<l0line*(1.0+zref)+700. and float(ents[0])>l0line*(1.0+zref)-700.:
            speclam.append(float(ents[0]))
            specval.append(float(ents[1]))
            specunc.append(float(ents[2]))
    infUL.close()
    if l0line*(1.0+zref)<min(speclam) or l0line*(1.0+zref)>max(speclam):
        return 0
    indclosest=-1
    dlclosest=100000.
    for p in range(len(speclam)):
        if math.fabs(speclam[p]-l0line*(1.0+zref))<dlclosest:
            indclosest=p
            dlclosest=math.fabs(speclam[p]-l0line*(1.0+zref))
    if indclosest==-1 or indclosest==0 or indclosest==len(speclam)-1:
        return 0
    else:
        dlam=speclam[indclosest]-speclam[indclosest-1]
    NpixFWHM=int(FWHMref*(1.0+zref)/(dlam))
    #print FWHMref,(1.0+zref),dlam
    if NpixFWHM%2==0:
        NpHWHM=NpixFWHM/2
    else:
        NpHWHM=(NpixFWHM-1)/2
    Npix=2*NpHWHM+1
    if indclosest-NpHWHM<0 or indclosest+NpHWHM>=len(speclam):
        print "Bad FWHM %s" % (specpath)
        return 0
    p=indclosest-NpHWHM
    sumsig2=0.0
    N=0
    while p<=indclosest+NpHWHM:
        #print p, len(specunc), NpHWHM
        sumsig2=sumsig2+specunc[p]**2
        N=N+1
        p=p+1
    #print N-Npix
    sigLINE=FWHMref*(1.0+zref)/(2.35*(dlam))
    nsig=(float(NpHWHM)+0.5)/sigLINE
    normalconst=quad(gauss,-1.0*nsig,nsig,full_output=0, epsabs=1.4899999999999999e-08, epsrel=1.4899999999999999e-08, limit=50, points=None, weight=None, wvar=None, wopts=None, maxp1=50, limlst=50)[0]
    #print normalconst
    return 3.0*sumsig2**0.5*dlam/normalconst
    
def go(linelistpath,outfilepath):
    fieldnos=[]
    objids=[]
    ras=[]
    decs=[]
    Js=[]
    Hs=[]
    redshifts=[]
    fluxes=[]
    fluxerrs=[]
    fwhms=[]
    fwhmerrs=[]
    lineids=[]

    infLL=open(linelistpath,'r')
    for line in infLL:
        if line[0]!='#':
            entries=line.split()
            fieldnos.append(entries[0])
            objids.append(entries[1])
            ras.append(entries[3])
            decs.append(entries[4])
            Js.append(entries[5])
            Hs.append(entries[6])
            redshifts.append(float(entries[10]))
            fluxes.append(float(entries[12]))
            fluxerrs.append(float(entries[13]))
            fwhms.append(float(entries[14]))
            fwhmerrs.append(float(entries[15]))
            lineids.append(entries[18])
    infLL.close()
    i=0
    outf=open(outfilepath,'w')
    print >>outf, "#1 Field ID\n#2 Object ID\n#3 RA\n#4 Dec\n#5 J mag\n#6 H mag\n#7 Redshift\n#8 [O II] flux\n#9 [O II] Flux Err"
    print >>outf, "#10 H beta Flux\n#11 H beta Flux Err\n#12 [O III] 5007 Flux\n#13 [O III] 5007 flux err\n#14 H alpha Flux\n#15 H alpha Flux Err"
    print >>outf, "#16 [S II] flux\n#17 [S II]flux err\n#18 [S III] 9069 flux\n#19 [S III] 9069 flux err\n#20 [S III] 9532 flux\n#21 [S III] 9532 flux err\n#22 He I 10830 flux\n#23 He I 10830 flux err"
    print >>outf, "#24 [O III] FWHM\n#25 H alpha FWHM\n#26 OII limiting flux\n#27 Hbeta limiting flux\n#28 OIII 5007 limiting flux\n#29 Halpha limiting flux\n#30 SII limiting flux\n#31 SIII 9069 limiting flux\n#32 SIII 9532 limiting flux"
    print >>outf, "#33 OII limiting flux 2\n#34 Hbeta limiting flux 2\n#35 OIII 5007 limiting flux 2\n#36 Halpha limiting flux 2\n#37 SII limiting flux 2\n#38 SIII 9069 limiting flux 2\n#39 SIII 9532 limiting flux 2"

    while i<len(fieldnos):
        if i==0 or fieldnos[i]!=fieldnos[i-1] or objids[i]!=objids[i-1]:
            o2=0.0
            o2err=99.0
            hb=0.0
            hberr=99.0
            o3=0.0
            o3err=99.0
            ha=0.0
            haerr=99.0
            s2=0.0
            s2err=99.0
            s31=0.0
            s31err=99.0
            s32=0.0
            s32err=99.0
            he=0.0
            heerr=99.0
            ra=ras[i]
            dec=decs[i]
            par=fieldnos[i]
            obj=objids[i]
            Jmag=Js[i]
            Hmag=Hs[i]
            z=redshifts[i]
            hafwhm=0.0
            o3fwhm=0.0
            hbfwhm=0.0
            o2fwhm=0.0
        if lineids[i]=='OII':
            o2=fluxes[i]
            o2err=fluxerrs[i]
            o2fwhm=fwhms[i]
        elif lineids[i]=='H_beta':
            hb=fluxes[i]
            hberr=fluxerrs[i]
            hbfwhm=fwhms[i]
        elif lineids[i]=='OIII_4959+5007':
            o3=0.75*fluxes[i]
            o3err=0.75*fluxerrs[i]
            o3fwhm=fwhms[i]
        elif lineids[i]=='H_alpha':
            ha=fluxes[i]
            haerr=fluxerrs[i]
            hafwhm=fwhms[i]
        elif lineids[i]=='SII':
            s2=fluxes[i]
            s2err=fluxerrs[i]
        elif lineids[i]=='SIII_9069':
            s31=fluxes[i]
            s31err=fluxerrs[i]
        elif lineids[i]=='SIII_9532':
            s32=fluxes[i]
            s32err=fluxerrs[i]
        elif lineids[i]=='HeI_10830':
            he=fluxes[i]
            heerr=fluxerrs[i]
        #Get Upper limits!
        if i==len(lineids)-1 or objids[i+1]!=objids[i]:
            limlams=[3727.0,4861.0,5007.0,6563.0,6724.0,9069.0,9532.0]
            lims=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            lims2=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            #lims=[]
            for q in range(len(limlams)):
                spath='../Spectra/Par'+str(par)+'_BEAM_'+str(obj)+'A.dat'
                limsL=0.0
                limsR=0.0
                limsL=getUL(limlams[q],limlams[q]-50.,z,max([hafwhm,o3fwhm]),spath)*3.0
                limsR=getUL(limlams[q],limlams[q]+50.,z,max([hafwhm,o3fwhm]),spath)*3.0
                if limsL!=0.0 and limsR!=0.0:
                    lims[q]=min([limsL,limsR])
                else:
                    lims[q]=max([limsL,limsR])
                path102='../Spectra/Par'+str(par)+'_G102_BEAM_'+str(obj)+'A.dat'
                path141='../Spectra/Par'+str(par)+'_G141_BEAM_'+str(obj)+'A.dat'
                if hafwhm>0.0 and o3fwhm>0.0:
                    FWHM=min([hafwhm,o3fwhm])
                else:
                    FWHM=max([hafwhm,o3fwhm])
                if limlams[q]*(1.0+z)>8000. and limlams[q]*(1.0+z)<11500. and os.path.exists(path102)==1:
                    lims2[q]=getULerrspec(limlams[q],z,FWHM,path102)
                elif limlams[q]*(1.0+z)>11500. and limlams[q]*(1.0+z)<17000. and os.path.exists(path141)==1:
                    lims2[q]=getULerrspec(limlams[q],z,FWHM,path141)
                #lims.append(getUL(limlams[q],limlams[q]-50.,z,max([hafwhm,o3fwhm]),spath))
                #if lims[q]>0.0:
                #    print lims[q]
            print >>outf, "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e" % (par,obj,ra,dec,Jmag,Hmag,z,o2,o2err,hb,hberr,o3,o3err,ha,haerr,s2,s2err,s31,s31err,s32,s32err,he,heerr,o3fwhm,hafwhm,lims[0],lims[1],lims[2],lims[3],lims[4],lims[5],lims[6],lims2[0],lims2[1],lims2[2],lims2[3],lims2[4],lims2[5],lims2[6])
        #Print to file!
        i=i+1
    outf.close()
