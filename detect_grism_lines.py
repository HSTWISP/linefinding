#################################################################################################
# detect_grism_lines.py v1.0
# Last Updated July 16, 2015 by Sophia Dai
# Last Update: added two lines in *previous* line 80 and 203
# to skip the lines starting with '#' int the spectra.dat files
#   if line.startswith('#'):
#            continue
#
# Written by Dr. Alaina Henry, UCSB and Nathaniel Ross, UCLA
#            ahenry@physics.ucsb.edu    nross@astro.ucla.edu
#
# to run this code, compile lists of spectra in the G102 and G141 in which you would like to
# automatically detect emission lines.  Note that you should not use the overlaid spectra 
# (i.e. both the g141 and g102 spectra together in a single file) as this would cause the code to
# fail and/or find spurious lines in the g141 in the region of overlap.  You should also create an 
# output directory where you would like the results of this code to be stored.
# The command line to run this routine is as follows:
#
# python detect_grism_lines.py [g102_list] g141_list output_directory output_line_list_file
#
# Note that if you pass only three arguments, this program assumes you want to fit G141 spectra only.
# If you pass it four arguments, it will fit both the G102 and G141 spectra.  Note also that the 
# output directory name should not have a trailing '/'.  Other possible sources of error include:
# having a blank line at the end of one of your lists of spectra, incorrectly named entries in the
# lists of spectra.
#################################################################################################
import os
import distutils
import numpy as np
import scipy
import matplotlib.pyplot as plt
#import matplotlib
#from numpy import r_, sin
from scipy.interpolate import spline
from distutils.sysconfig import *
import math

from pyraf import iraf # May need link to login.cl file in the working directory
#from iraf import noao

iraf.onedspec()
iraf.unlearn('scombine')
iraf.unlearn('continuum')
iraf.unlearn('rspectext')

cont_function='spline3'
cont_order='3'


###############################################################################
# This function will flag and remove from the queue those spectra which contain
# too many bad pixels to properly fit the continuum.
def reject_bad_spectra (spec_list,minlam,maxlam,outfile):
    good_spec=[]
    bad_spec=[]
    for infile in spec_list:
        l=len(infile)
        name=infile[0:l-1]
#        name = name[0:l-6]+'.dat'
        inspec = open(name)
        nbadpix = 0
        for line in inspec:
            entries = line.split()
            if float(entries[1])==0.00 and float(entries[0])>=minlam and float(entries[0])<=maxlam:
                nbadpix=nbadpix+1
        inspec.close()
        if nbadpix>33:
            bad_spec.append(infile)
        else:
            good_spec.append(infile)
    if len(bad_spec)>0:
        outf=open(outfile,'w')
        for badf in bad_spec:
            outf.write(badf)
        outf.close
    return good_spec
##############################################################################################################
def get_spec_limits (specdatfile,grismmin,grismmax):
    specfile=open(specdatfile,'r')
    speclam=[]
    speccut=[]
    minlam=grismmin
    maxlam=grismmax
    for line in specfile:
        if line.startswith('#'):
            continue
        entries = line.split()
        wavelen = float(entries[0])
        speclam.append(wavelen)
        if len(entries)>=5:
            cutoff=float(entries[4])
            if entries[2]=='nan' or entries[2]=='NaN' or entries[2]=='NAN' or entries[2]=='-nan' or entries[2]=='-NaN' or entries[2]=='-NAN':
                cutoff=2.0
            speccut.append(cutoff)
    specfile.close()
    if len(speccut)==0:
        result=[grismmin,grismmax]
        return result
    
    nspec=len(speclam)
    midpt=nspec/2
    minset=0
    for t in range(nspec):
#        if t<midpt and speclam[t]>=minlam and speccut[t]==2.0:
        if minset==0 and speclam[t]>=minlam and speccut[t]==2.0 and t<(nspec-4):
            minlam=speclam[t+3]
        if minset==0 and speclam[t]>=minlam and speccut[t]==0.0:
            minset=1
#        if t>midpt and speclam[t]<=maxlam and speccut[t]==2.0:
        if minset==1 and speclam[t]<=maxlam and speccut[t]==2.0:
            maxlam=speclam[t-3]
    if minlam>grismmax:
        result=[0.0,0.0]
    elif minlam<maxlam:
        result=[minlam,maxlam]
    else:
        result=[grismmin,grismmax]
    return result

###############################################################################
# This function will fit the continuum to a spectrum and output the result in 
# both a .fits and a .dat format
def fit_cont (specdat,outputdir,minlam,maxlam):
    p=len(specdat)
    #specdat=specdat[0:p-1]
    specfits = specdat[0:p-4] + '.fits'
    outfile_clip = outputdir +  '/' + specdat[0:p-4] + '.clip.fits'
    outfile_cont_fit_fits = outputdir + '/'  + specdat[0:p-4] + '.contfit.fits'
    outfile_cont_fit_dat = outputdir + '/'  + specdat[0:p-4] + '.contfit.dat'
    if os.path.exists(specfits) ==1:
        os.unlink(specfits)
    if os.path.exists(outfile_clip) ==1:
        os.unlink(outfile_clip)
    if os.path.exists(outfile_cont_fit_fits) ==1:
        os.unlink(outfile_cont_fit_fits)
    if os.path.exists(outfile_cont_fit_dat) ==1:
        os.unlink(outfile_cont_fit_dat)

    iraf.rspectext(specdat,specfits,dtype='nonlinear',crval='INDEF',cdelt1='INDEF')
#    print 'Steve Holt! 1'
    iraf.scombine(specfits, outfile_clip, w1=minlam, w2=maxlam)
#    print 'Steve Holt! 2'
    iraf.continuum(outfile_clip, outfile_cont_fit_fits, type = 'fit', interactive='no', function=cont_function, order=cont_order, low_reject='2', high_reject='2', naverage='4') # Fit continuum
#    print 'Steve Holt! 3'
    iraf.wspectext(outfile_cont_fit_fits, outfile_cont_fit_dat, header = 'no') # Output in .dat format
#    print 'Steve Holt! 4'
    if os.path.exists(specfits) ==1:
        os.unlink(specfits)               #Added 08/22/2012 to delete extraneous data outputs--NRR
    if os.path.exists(outfile_clip) ==1:
        os.unlink(outfile_clip)           #Added 08/22/2012 to delete extraneous data outputs--NRR
    if os.path.exists(outfile_cont_fit_fits) ==1:
        os.unlink(outfile_cont_fit_fits)  #Added 08/22/2012 to delete extraneous data outputs--NRR

# End Function fit_cont
###############################################################################

###############################################################################
# This function will read in the spectrum and continuum fit, subtract them, 
# measure the s/n, search for significant lines, and return their location.
# Also plots the spectrum, continuum fit, continuum-subtracted spectrum.
def find_lines (specdat,outputdir,minlam,maxlam,outlist,chi2outlist):
    minew=10.0   # Set minimum EW for acceptable line candidate
    maxcontig=20 # Set maximum pixel width for acceptable line candidate.
    p=len(specdat)
    lambda_min=minlam
    lambda_max=maxlam
    #specdat=specdat[0:p-1]
    beamname=specdat.split('_')
    fieldid=beamname[0]
    fieldid=fieldid[3:]
    beamstr=beamname[3]
    filtid=beamname[1]
    l1=len(beamstr)
    beamno=beamstr[0:l1-5]
    infile_continuum = outputdir + '/'  + specdat[0:p-4] + '.contfit.dat'
    infile_spectrum = specdat
    outfile_cont_fit_eps = outputdir + '/'  + specdat[0:p-4] + '_contfit.eps'
    if os.path.exists(outfile_cont_fit_eps) ==1:
        os.unlink(outfile_cont_fit_eps)

    continuum_lam=[]
    continuum_val=[]
    

	# Read in continuum values
#    print 'Steve Holt! 5'
    continuumin=open(infile_continuum,'r')
    for line in continuumin:
	entries=line.split()
	if float(entries[0])>=lambda_min and float(entries[0])<=lambda_max:
		continuum_lam.append(float(entries[0]))
		continuum_val.append(float(entries[1]))
		#print continuum_lam[j],continuum_val[j] #Check Read-in
	#		j=j+1
    continuumin.close()
	
	#j=0
    spectrum_lam=[]
    spectrum_val=[]
    spectrum_unc=[]
    spectrum_con=[]
    spectrum_zer=[]
	#Read in spectrum values
#    print 'Steve Holt! 6'
    spectrumin=open(infile_spectrum,'r')
    #bad_pix=0
    for line in spectrumin:
        if line.startswith('#'):
            continue
        entries=line.split()
        if float(entries[0])>=lambda_min and float(entries[0])<=lambda_max:
            spectrum_lam.append(float(entries[0]))
            spectrum_val.append(float(entries[1]))
            if entries[2]=='nan' or entries[2]=='NAN' or entries[2]=='NaN' or entries[2]=='-NaN':
                spectrum_unc.append(0.00)
                #bad_pix=bad_pix+1 # Flag bad pixels for possible rejection of problem spectra
            else:
                spectrum_unc.append(float(entries[2]))
	#		j=j+1
            if len(entries)>=4:
                spectrum_con.append(float(entries[3]))
            if len(entries)>=5:
                spectrum_zer.append(float(entries[4]))
    spectrumin.close()
	#print bad_pix
#    print 'Steve Holt! 6.1'
    cont_lam = np.array(continuum_lam)
    cont_val = np.array(continuum_val)
    spec_lam = np.array(spectrum_lam)
    spec_val = np.array(spectrum_val)
    spec_unc = np.array(spectrum_unc)
    if len(spectrum_con)>0:
        spec_con = np.array(spectrum_con)
    if len(spectrum_zer)>0:
        spec_zer = np.array(spectrum_zer)
#    print 'Steve Holt! 6.3'
    print len(cont_lam), len(cont_val), len(spec_lam)
    cont_interp_val = scipy.interpolate.spline(cont_lam,cont_val,spec_lam)
#    print 'Steve Holt! 6.5'
    cont_subtracted_spec = spec_val - cont_interp_val
    #noise = np.std(cont_subtracted_spec, axis=None, dtype=None, out=None, ddof=0)
#    print 'Steve Holt! 7'
    chi2=0.00
    nbad=0.00
    for k in range(len(cont_subtracted_spec)):
        if spec_unc[k]>0.00:
            chi2=chi2+(cont_subtracted_spec[k]**2/spec_unc[k]**2)
        else:
            nbad=nbad+1.00
    
    if cont_function=='spline3':
        n_params=4.0*float(cont_order) - 3.0*(float(cont_order)-1.0)
    elif cont_function=='spline1':
        n_params=2.0*float(cont_order) - (float(cont_order)-1.0)
    else:
        n_params=float(cont_order)+1.0
    
    red_chi2=chi2/(float(len(cont_subtracted_spec))-nbad-n_params)
    
    #median_unc=np.median(spec_unc)
#    print 'Steve Holt! 8'
    s_to_n=[]
    n1=len(spec_lam)
    n2=len(cont_lam)
    for j in range(n1):
        if spectrum_unc[j]>0.0:
            #s_to_n.append((spectrum_val[j]-interp_val[j])/spectrum_unc[j])
            #s_to_n.append(cont_subtracted_spec[j]/(noise**2.0+(0.01*spec_unc[j])**2)**0.5)
            s_to_n.append(cont_subtracted_spec[j]/spec_unc[j])
            #s_to_n.append(cont_subtracted_spec[j]/(spec_unc[j]*noise/median_unc))
            #s_to_n.append(cont_subtracted_spec[j]/(spec_unc[j]*(red_chi2)**0.5))
        else:
            s_to_n.append(0.00)
    flags=0
    if red_chi2>20.0:
        flags=flags+1
    print >>chi2outlist,"%s\t%f\t%f" % (specdat,red_chi2,n_params)
#    print 'Steve Holt! 9'
    ########################################################################
#    n_sig=2.88 #Number of sigma above which a pixel is considered significant
    n_sig=1.73
    n_sig_2pix=3.54 #For a 2-pixel line, ensure that > 5-sigma significant
    n_thresh=3 #Number of contiguous significant pixels to define a line
    ########################################################################
    n_contig=0
    #flux=0
    #contEW=[]
    index=0
    print infile_spectrum
    #Find strings of contiguous pixels longer than n_thresh with s/n greater than n_sig
    linelam=[]
    lineval=[]
    for j in range(n1):
       	if s_to_n[j]>=n_sig:
            n_contig=n_contig+1 #Start counting!			
        if s_to_n[j]<n_sig and n_contig>=n_thresh:
            index=j - n_contig/2 - 1
            #print "Significant line of %i contiguous pixels at %f Angs." % (n_contig,spec_lam[index]) #Found a line!
            tot_noise=0.00
            tot_signal=0.00
            ewest=0.0
            ewcont=[]
            dlam=spec_lam[j-n_contig+1]-spec_lam[j-n_contig]
            tot_ston=0.00
            zero_order=0
            cutoff=0
            tot_contam=0.00
            tot_spect=0.00
            #print spec_lam[index]
            for i in range(n_contig):
                
                tot_signal=tot_signal+cont_subtracted_spec[j-n_contig+i]
                #print tot_signal,cont_subtracted_spec[j-n_contig+i]
                ewcont.append(cont_interp_val[j-n_contig+i])
                tot_noise = tot_noise + (spec_unc[j-n_contig+i]**2.0)*red_chi2
                if len(spectrum_con)>0:
                    tot_contam = tot_contam + spectrum_con[j-n_contig+i]
                    tot_spect = tot_spect + spec_val[j-n_contig+i]
                if len(spectrum_zer)>0 and spectrum_zer[j-n_contig+i]==1.0:
                    zero_order=1
            if (j-n_contig)<10:
                cutrange=j+10
                ncut=j-n_contig
            elif (n1-j)<10:
                cutrange=n1-j+10+n_contig
                ncut=10
            else:
                cutrange=20+n_contig
                ncut=10
            for k in range(cutrange):
                if len(spectrum_zer)>0 and spectrum_zer[j-n_contig+k-ncut-1]==2.0:
                    cutoff=1
            tot_noise=tot_noise**0.5
            tot_ston=tot_signal/tot_noise
            #print spec_lam[index],tot_signal,tot_noise,tot_ston
            ewest=tot_signal*dlam/np.median(np.array(ewcont))
            if n_contig>20:
                flags=flags+2
            if zero_order==1:
                flags=flags+4
            if cutoff==1:
                flags=flags+8
            if len(spectrum_con)>0:
                contamination_level=tot_contam/tot_spect
                if contamination_level>=0.20:
                    flags=flags+16
            #if n_contig<maxcontig and ewest>minew: # Check if line too large or too low EW to be real
            if math.fabs(ewest)>minew: # Check if line too low EW to be real
                print >>outlist, "%s\t%s\t%s\t%f\t%i\t%f\t%i " % (fieldid,filtid,beamno,spec_lam[index],n_contig,tot_ston,flags)
                #print >>outlist, "%s\t%s\t%s\t%f\t%i\t%f\t%i\tEW%.3f\tSIG%.3e\tDL%.3e\tCONT%.3e " % (fieldid,filtid,beamno,spec_lam[index],n_contig,tot_ston,flags,ewest,tot_signal,dlam,np.median(np.array(ewcont)))
            linelam.append(spec_lam[index])
            lineval.append(spec_val[index])
## Eliminated 2-pixel criterion for G102 on 2011/09/19 per James Colbert's request
        if s_to_n[j]<n_sig and n_contig==2 and (s_to_n[j-2]>=n_sig_2pix and s_to_n[j-1]>=n_sig_2pix) and spec_lam[j]>11500.:
            tot_signal=cont_subtracted_spec[j-2]+cont_subtracted_spec[j-1]
            ewcont=[cont_interp_val[j-2],cont_interp_val[j-1]]
            tot_noise=((spec_unc[j-2]**2.0+spec_unc[j-1]**2.0)*red_chi2)**0.5
            tot_ston=tot_signal/tot_noise
            dlam=spec_lam[j-1]-spec_lam[j-2]
            ewest=tot_signal*dlam/np.median(np.array(ewcont))
            if math.fabs(ewest)>minew:
                print >>outlist, "%s\t%s\t%s\t%f\t%i\t%f\t%i " % (fieldid,filtid,beamno,spec_lam[j-2],n_contig,tot_ston,flags)
                #print >>outlist, "%s\t%s\t%s\t%f\t%i\t%f\t%i\tEW%.3f\tSIG%.3e\tDL%.3e\tCONT%.3e " % (fieldid,filtid,beamno,spec_lam[index],n_contig,tot_ston,flags,ewest,tot_signal,dlam,np.median(np.array(ewcont)))
            linelam.append(spec_lam[j-2])
            lineval.append(spec_val[j-2])
        if s_to_n[j]<n_sig:
            n_contig=0 #Stop counting!

#    print 'Steve Holt! 10'
    ###############################################################################
    #Plot Spectra:
    #lam = np.array(spectrum_lam)
    #spec = np.array(spectrum_val)
    #unc = np.array(spectrum_unc)
    #fit_lam = np.array(continuum_lam)
    #it_spec = np.array(continuum_val)
	
    #crude automatic y scaling.
    #ymin = spec_val.min()
    #ymin = cont_subtracted_spec.min()
    #ymax = 1.5*spec_val.max()
	
    #title = specdat[0:p-5] 

    #clight=3.0*(10.0**8.0)*(10.0**10.0) # Speed of light in units of Angs*Hz

    #l_lam = np.array(linelam)
    #l_val = np.array(lineval)
    #l_val = l_val*(l_lam**2)/clight
    
    #spec_fnu = spec_val*(spec_lam**2)/clight
    #cont_fnu = cont_val*(cont_lam**2)/clight
    #contam_fnu = spec_con*(spec_lam**2)/clight
    #zero_ord = spec_zer * spec_fnu

    #ymin=spec_fnu.min()
    #ymax=1.5*spec_fnu.max()
    
    ##### The following was commented out on 08/22/2012 to remove extraneous data output--NRR
    #plt.clf()
    #plt.ioff()
    #plt.figure(1,figsize=(7, 10))
    #plt.subplot(211)
    ##plt.plot(spec_lam, spec_val, 'k', cont_lam, cont_val, 'r', spec_lam, cont_subtracted_spec,'b',spec_lam,spec_unc,'g')
    ##specplt = plt.plot(spec_lam, spec_val, 'k')
    #specplt = plt.plot(spec_lam, spec_fnu, 'k',spec_lam,contam_fnu,'r',spec_lam,zero_ord,'m')
    #plt.setp(specplt, ls = 'steps')
    ##curveplt = plt.plot(spec_lam, cont_subtracted_spec, 'b',  cont_lam, cont_val, 'r')
    ##curveplt = plt.plot(cont_lam, cont_val, 'r')
    #curveplt = plt.plot(cont_lam, cont_fnu, 'b')
    #plt.xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
    #plt.ylabel(r'F$_\nu$ ergs s$^{-1}$ cm$^{-2}$ Hz$^{-1}$', size='xx-large')
    #plt.ylim([ymin, ymax])
    #if len(l_lam)>0:
    #    linepts=plt.scatter(l_lam,l_val,s=80,c='g',marker='o', cmap=None, norm=None,vmin=None, vmax=None, alpha=1, linewidths=None,verts=None)
    #plt.title(title)
    ###plt.show()  # display plot on your screen!
    ###plt.savefig(outfile_cont_fit_eps)
    #
    #plt.subplot(212)
    #threshline_x=np.array([minlam,maxlam])
    #threshline_y=np.array([n_sig,n_sig])
    #s_to_n_eps = outputdir + '/' + specdat[0:p-5] + '_signal_to_noise.eps'
    ##plt.clf()
    ##plt.plot(spec_lam, s_to_n, 'k')
    #stonspec= plt.plot(spec_lam,s_to_n,'k',threshline_x,threshline_y,'r')
    #plt.setp(stonspec, ls = 'steps')
    #plt.xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
    #plt.ylabel(r'S/N')
    ##plt.title(title)
    ###plt.savefig(s_to_n_eps)
    #plt.savefig(outfile_cont_fit_eps) 
#    print 'Steve Holt! 11'
    ##############################################################################
# End line detection function find_lines()
####################################################################################


def get_unique(somelist):
    somelist.sort()
    resultlist=[]
    for jk in range(len(somelist)):
        if jk==0 or somelist[jk]!=somelist[jk-1]:
            resultlist.append(somelist[jk])
    return resultlist

###############################################################################

def merge_emitters(linelistname):
    llfile=linelistname
    parts=llfile.split('.')
    pdfname=parts[0]+'.pdf'
    all_objs=[]
    llin=open(linelistname,'r')
    for line in llin:
        entries=line.split()
        all_objs.append(entries[2])
        parno=entries[0]
    llin.close()
    uniq_objs=get_unique(all_objs)
    mergestr=''
    for obj in uniq_objs:
        g102e='linelist/Par' + str(parno) + '_G102_BEAM_'+ str(obj) + 'A_contfit.eps'
        g141e='linelist/Par' + str(parno) + '_G141_BEAM_'+ str(obj) + 'A_contfit.eps'
        if os.path.exists(g102e)==1:
            mergestr=mergestr + ' ' + g102e
        if os.path.exists(g141e)==1:
            mergestr=mergestr + ' ' + g141e
    cmd = 'gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=' + pdfname + mergestr
    os.system(cmd)


###############################################################################
# Begin body of the program:
###############################################################################

#define the wavelength clipping, to avoid the ends of the spectra
# we might need to change these. 
lambda_min_g102 = 8000.
lambda_max_g102 = 1.15e4
lambda_min_g141 = 1.1e4
lambda_max_g141 = 1.67e4

######################################################################

outdir='linelist'

if os.path.exists('all_spec.list')==1:
    os.unlink('all_spec.list')

os.system('ls *G102_BEAM_*A.dat > all_spec.list')
os.system('ls *G141_BEAM_*A.dat >> all_spec.list')

existsG102=0
existsG141=0

allfiles=[]
allbeamnos=[]

readspecfiles=open('all_spec.list','r')
for fname in readspecfiles:
    
    filename=fname
    flen=len(fname)
    filename=filename[0:flen-1]
    allfiles.append(filename)
    fparts=filename.split('_')
    if existsG102==0 and fparts[1]=='G102':
        existsG102=1
        print 'G102 files detected'
    if existsG141==0 and fparts[1]=='G141':
        existsG141=1
        print 'G141 files detected'
    beamlen=len(fparts[3])
    beamno=fparts[3]
    allbeamnos.append(int(beamno[0:beamlen-5]))
readspecfiles.close()

if os.path.exists('linelist/')==0:
    os.system('mkdir linelist')
else:
    os.system('rm -f linelist/*.eps')

prototype=allfiles[0].split('_')
parid=prototype[0]
linelistfile='linelist/' + parid + 'lines.dat'
chi2outname='linelist/' + parid + '_reducedchi2.dat'
linelistout=open(linelistfile,'w')
chi2out=open(chi2outname,'w')

print "I detect %i files" % (len(allbeamnos))

uniquelist=get_unique(allbeamnos)
print "I detect %i unique beam id numbers" % (len(uniquelist))

for beam in uniquelist:
    beam102=prototype[0] + '_G102_BEAM_' + str(beam) + 'A.dat'
    beam141=prototype[0] + '_G141_BEAM_' + str(beam) + 'A.dat'
    if os.path.exists(beam102)==1:
        #print beam102
        limitsg102=get_spec_limits(beam102,lambda_min_g102,lambda_max_g102)
        if limitsg102[0]!=0.0 and limitsg102[1]!=0.0 and (limitsg102[1]-limitsg102[0])>1000.0:
#        fit_cont(beam102,outdir,lambda_min_g102,lambda_max_g102)
            fit_cont(beam102,outdir,limitsg102[0],limitsg102[1])
#        find_lines(beam102,outdir,lambda_min_g102,lambda_max_g102,linelistout,chi2out)
            find_lines(beam102,outdir,limitsg102[0],limitsg102[1],linelistout,chi2out)
    if os.path.exists(beam141)==1:
        limitsg141=get_spec_limits(beam141,lambda_min_g141,lambda_max_g141)
#        fit_cont(beam141,outdir,lambda_min_g141,lambda_max_g141)
        if limitsg141[0]!=0.0 and limitsg141[1]!=0.0 and (limitsg141[1]-limitsg141[0])>1000.0:
            fit_cont(beam141,outdir,limitsg141[0],limitsg141[1])
#        find_lines(beam141,outdir,lambda_min_g141,lambda_max_g141,linelistout,chi2out)
            find_lines(beam141,outdir,limitsg141[0],limitsg141[1],linelistout,chi2out)
       

    
linelistout.close()
chi2out.close()
#if os.path.exists('merged.pdf')==1:
#    os.unlink('merged.pdf')
#print 'Merging pdfs...this may take awhile'
#cmd = 'gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged.pdf ' + outdir + '/*.eps'
#os.system(cmd)
##merge_emitters(linelistfile) ## Removed on 08/22/2012 to remove extraneous data output--NRR
print 'Complete!'

