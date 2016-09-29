pro matchextraction, file1, file2, pref3, pref4, file5, prefix

	snr_limit=3.0
        snr_limit2=2.0

        file3a=pref3+'.dat'
        file4a=pref4+'.dat'
        file3=pref3+'_contam.dat'
        file4=pref4+'_contam.dat'

	readcol,file1,par,obj,cont,ra,dec,jmag,hmag,major,minor,rot,redshift,redshiftE,flux,fluxE,fwhm,fwhmE,EW,EWE,ID,$
                format='i,i,i,d,d,f,f,f,f,f,f,f,f,f,f,f,f,f,a',skipline=19
	readcol,file2,par2,obj2,cont2,ra2,dec2,jmag2,hmag2,major2,minor2,rot2,redshift2,redshiftE2,flux2,fluxE2,fwhm2,fwhmE2,EW2,EWE2,ID2,$
                format='i,i,i,d,d,f,f,f,f,f,f,f,f,f,f,f,f,f,a',skipline=19
  

        readcol,file3,par3,format='i'
        readcol,file4,par4,format='i'
        ltt3=n_elements(par3)
        ltt4=n_elements(par4)
        readcol,file3a,par3a,format='i'
        readcol,file4a,par4a,format='i'
        lta3=n_elements(par3a)
        lta4=n_elements(par4a)

        grism='GGGG'
        stuff=0
        size=0.0
        flag=1
        contam=0
        z=0.0
        comment='ggggggggggggggggggggggggggggggggg'

        par3a=intarr(lta3)
        obj3a=intarr(lta3)
        comment3a=strarr(lta3)
        readcol,file3a,testpar,junk,testobj,testwave,teststuff,testsize,testflag,format='i,a,i,f,i,f,i'
        if testpar(0) lt 10 then begin
           parnum=2
        endif else if testpar(0) lt 100 then begin
           parnum=3
        endif else parnum=4
        insertpar=strcompress(string(parnum),/REMOVE_ALL)
        openr,1,file3a
        for i=0,lta3-1 do begin
           danum=5
           if testobj(i) lt 10 then begin
              danum=2
           endif else if testobj(i) lt 100 then begin
              danum=3
           endif else if testobj(i) lt 1000 then begin
              danum=4
           endif
           danum2=3
           if teststuff(i) lt 10 then danum2=2
           danum3=2
           if testflag(i) ge 10 then danum3=3
           danum4=12
           if testwave(i) ge 10000 then danum4=13
           insert1=strcompress(string(danum),/REMOVE_ALL)
           insert2=strcompress(string(danum2),/REMOVE_ALL)
           insert3=strcompress(string(danum3),/REMOVE_ALL)
           insert4=strcompress(string(danum4),/REMOVE_ALL)
           readf,1,dapar,grism,daobj,wave,stuff,size,flag,contam,z,comment,format='(I'+insertpar+',A5,I'+insert1+',F'+insert4+'.5,I'+insert2+',F9.6,I'+insert3+',I2,F9.6,A100)'
           test=strcompress(comment,/remove_all)
          ; print,daobj,comment
           if test eq '0' OR test eq '1' OR test eq '2' OR test eq '3' OR test eq '4' OR test eq '5' OR test eq '6' $
              OR test eq '7' OR test eq '8' OR test eq '9' OR test eq '00' then comment=' '
           par3a(i)=dapar
           obj3a(i)=daobj
           comment3a(i)=strcompress(comment)
           ;print,comment3a(i),obj3a(i)
        endfor
        close,1     
        par4a=intarr(lta4)
        obj4a=intarr(lta4)
        comment4a=strarr(lta4)
        readcol,file4a,testpar,junk,testobj,testwave,teststuff,testsize,testflag,format='i,a,i,f,i,f,i'
        openr,1,file4a
        for i=0,lta4-1 do begin
           danum=5
           if testobj(i) lt 10 then begin
              danum=2
           endif else if testobj(i) lt 100 then begin
              danum=3
           endif else if testobj(i) lt 1000 then begin
              danum=4
           endif
           danum2=3
           if teststuff(i) lt 10 then danum2=2
           danum3=2
           if testflag(i) ge 10 then danum3=3
           danum4=12
           if testwave(i) ge 10000 then danum4=13
           insert1=strcompress(string(danum),/REMOVE_ALL)
           insert2=strcompress(string(danum2),/REMOVE_ALL)
           insert3=strcompress(string(danum3),/REMOVE_ALL)
           insert4=strcompress(string(danum4),/REMOVE_ALL)
           readf,1,dapar,grism,daobj,wave,stuff,size,flag,contam,z,comment,format='(I'+insertpar+',A5,I'+insert1+',F'+insert4+'.5,I'+insert2+',F9.6,I'+insert3+',I2,F9.6,A100)'
           test=strcompress(comment,/remove_all)
           if test eq '0' OR test eq '1' OR test eq '2' OR test eq '3' OR test eq '4' OR test eq '5' OR test eq '6' $
              OR test eq '7' OR test eq '8' OR test eq '9' OR test eq '00' then comment=' '
           par4a(i)=dapar
           obj4a(i)=daobj
           comment4a(i)=strcompress(comment)
        endfor
        close,1     
        

        if ltt3 gt 0 then begin
        par3=intarr(ltt3)
        obj3=intarr(ltt3)
        comment3=strarr(ltt3)
        readcol,file3,testpar,junk,testobj,testwave,teststuff,testsize,testflag,format='i,a,i,f,i,f,i'
        openr,1,file3
        for i=0,ltt3-1 do begin
           danum=5
           if testobj(i) lt 10 then begin
              danum=2
           endif else if testobj(i) lt 100 then begin
              danum=3
           endif else if testobj(i) lt 1000 then begin
              danum=4
           endif
           danum2=3
           if teststuff(i) lt 10 then danum2=2
           danum3=2
           if testflag(i) ge 10 then danum3=3
           danum4=12
           if testwave(i) ge 10000 then danum4=13
           insert1=strcompress(string(danum),/REMOVE_ALL)
           insert2=strcompress(string(danum2),/REMOVE_ALL)
           insert3=strcompress(string(danum3),/REMOVE_ALL)
           insert4=strcompress(string(danum4),/REMOVE_ALL)
           readf,1,dapar,grism,daobj,wave,stuff,size,flag,contam,z,comment,format='(I'+insertpar+',A5,I'+insert1+',F'+insert4+'.5,I'+insert2+',F9.6,I'+insert3+',I2,F9.6,A100)'
           test=strcompress(comment,/remove_all)
           if test eq '0' OR test eq '00' then comment=' '
           par3(i)=dapar
           obj3(i)=daobj
           comment3(i)=strcompress(comment)
        endfor
        close,1      
        endif

        if ltt4 gt 0 then begin
        par4=intarr(ltt4)
        obj4=intarr(ltt4)
        comment4=strarr(ltt4)
        readcol,file4,testpar,junk,testobj,testwave,teststuff,testsize,testflag,format='i,a,i,f,i,f,i'
        openr,1,file4
        for i=0,ltt4-1 do begin
           danum=5
           if testobj(i) lt 10 then begin
              danum=2
           endif else if testobj(i) lt 100 then begin
              danum=3
           endif else if testobj(i) lt 1000 then begin
              danum=4
           endif
           danum2=3
           if teststuff(i) lt 10 then danum2=2
           danum3=2
           if testflag(i) ge 10 then danum3=3
           danum4=12
           if testwave(i) ge 10000 then danum4=13
           insert1=strcompress(string(danum),/REMOVE_ALL)
           insert2=strcompress(string(danum2),/REMOVE_ALL)
           insert3=strcompress(string(danum3),/REMOVE_ALL)
           insert4=strcompress(string(danum4),/REMOVE_ALL)
           readf,1,dapar,grism,daobj,wave,stuff,size,flag,contam,z,comment,format='(I'+insertpar+',A5,I'+insert1+',F'+insert4+'.5,I'+insert2+',F9.6,I'+insert3+',I2,F9.6,A100)'
           test=strcompress(comment,/remove_all)
           if test eq '0' OR test eq '00' then comment=' '
           par4(i)=dapar
           obj4(i)=daobj
           comment4(i)=strcompress(comment)
        endfor
        close,1    
        endif

	lth=n_elements(par)
	lth2=n_elements(par2)
	good1=where(flux/fluxE gt snr_limit)
	good2=where(flux2/fluxE2 gt snr_limit)
	ok1=where(flux/fluxE gt snr_limit2)
	ok2=where(flux2/fluxE2 gt snr_limit2)
	print,n_elements(good1),n_elements(ok1),n_elements(good2),n_elements(ok2)

	gap='  '
	
; Make List With Only Brightest Lines For File1	
	match=intarr(lth)
	object=obj(0)
	snr=flux(0)/fluxE(0)
	theflux=flux(0)
	thefluxE=fluxE(0)
	theID=ID(0)
	thez=redshift(0)
	theEW=EW(0)
        thecont=cont(0)
	check=0
       
	openw,1,'linesort1.temp'
        if lth gt 1 then begin
	for i=1,lth-1 do begin
		if obj(i) eq object then begin
			if check eq 0 then begin
				match(i)=1.0
				check=1
			endif
			if flux(i)/fluxE(i) gt snr then begin
				theflux=flux(i)
				thefluxE=fluxE(i)
				snr=flux(i)/fluxE(i)
				theID=ID(i)
				thez=redshift(i)
        			theEW=EW(i)
                                thecont=cont(i)
			endif
		endif else begin
			if snr gt snr_limit then begin
                           printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
                           ;print,'GO',object,theflux,snr
                        endif
			check=0
			object=obj(i)
			theflux=flux(i) 
                        thefluxE=fluxE(i)
                        snr=flux(i)/fluxE(i)
		 	theID=ID(i)
			thez=redshift(i)
                        theEW=EW(i)
                        thecont=cont(i)
		endelse
	endfor
       if snr gt snr_limit then printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
       endif else begin
           if snr gt snr_limit then begin
              printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
           endif else print,'NO USABLE LINES IN LIST1!'
       endelse
       close,1
        

; Make List With Only Brightest Lines For File2 
        match2=intarr(lth2)
        object=obj2(0)
        snr=flux2(0)/fluxE2(0)
        theflux=flux2(0)
        thefluxE=fluxE2(0)
        theID=ID2(0)
        thez=redshift2(0)
        theEW=EW2(0)
        thecont=cont2(0)
        check=0

        openw,1,'linesort2.temp'
        if lth2 gt 1 then begin
        for i=1,lth2-1 do begin
                if obj2(i) eq object then begin
                        if check eq 0 then begin
                                match2(i)=1.0
                                check=1
                        endif
                        if flux2(i)/fluxE2(i) gt snr then begin
                                theflux=flux2(i)
                                thefluxE=fluxE2(i)
                                snr=flux2(i)/fluxE2(i)
                                theID=ID2(i)
                                thez=redshift2(i)
                                theEW=EW2(i)
                                thecont=cont2(i)
                        endif
                endif else begin
			if snr gt snr_limit then begin
                           printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
                          ; print,'GO',object,theflux,snr
                        endif
                        check=0
                        object=obj2(i)
                        theflux=flux2(i)
                        thefluxE=fluxE2(i)
                        snr=flux2(i)/fluxE2(i)
                        theID=ID2(i)
                        thez=redshift2(i)
                        theEW=EW2(i)
                        thecont=cont2(i)
                endelse
        endfor
        if snr gt snr_limit then printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
        endif else begin
           if snr gt snr_limit then begin
              printf,1,object,theflux,snr,theEW,thez,gap,thecont,gap,theID,FORMAT='(I6,E16.6,3F13.6,A6,I3,A6,A16)'
           endif else print,'NO USABLE LINES IN LIST2!
        endelse
	close,1

;Read temp files with brightest line only objects back in

	readcol,'linesort1.temp',object,theflux,snr,theEW,thez,thecont,theID,format='i,f,f,f,f,i,a'
	readcol,'linesort2.temp',object2,theflux2,snr2,theEW2,thez2,thecont2,theID2,format='i,f,f,f,f,i,a'

;Match the objects
	one=1
	two=2
        three=3
	zero=0
	na=-99.0

;These are the criteria for failed matches, presently 2% dif in redshift and 3 sigma difference in flux.
	delta=0.02
	delta2=3

	match,object,object2,a,b
	zratio=thez(a)/thez2(b)
	fluxratio=theflux(a)/theflux2(b)
        badID=where((zratio lt 1.0+delta AND zratio gt 1.0-delta) AND (theID(a) ne theID2(b)))
        badIDlth=n_elements(badID)
        nbadflux=intarr(n_elements(a))
        nbadlabel=strarr(n_elements(a))
        badfluxes=fltarr(n_elements(a),2)
        abadflux=intarr(n_elements(object))
        nbadflux(*)=0
        abadflux(*)=0
        if badIDlth(0) ne -1 then begin
        for i=0,badIDlth-1 do begin
           for j=0,lth2-1 do begin
              if object(a(badID(i))) eq obj2(j) AND theID(a(badID(i))) eq ID2(j) then begin
                 nfluxratio=theflux(a(badID(i)))/flux2(j)
                 print,'Check: ',object(a(badID(i))),obj2(j),theID(a(badID(i))),ID2(j),theflux(a(badID(i)))/flux2(j)
                 if nfluxratio lt (1.0+delta2*(1/snr(a(badID(i)))))^(-1) OR nfluxratio gt (1.0-delta2*(1/snr(a(badID(i)))))^(-1) then begin
                    print,'Found extra bad flux.'
                    nbadflux(badID(i))=1
                    nbadlabel(badID(i))=ID2(j)
                    badfluxes(badID(i),0)=theflux(a(badID(i)))
                    badfluxes(badID(i),1)=flux2(j)
                    abadflux(a(badID(i)))=1
                 endif
              endif
           endfor
           for j=0,lth-1 do begin
              if object2(b(badID(i))) eq obj(j) AND theID2(b(badID(i))) eq ID(j) then begin
                 nfluxratio=theflux2(b(badID(i)))/flux(j)
                 print,'Check: ',object2(b(badID(i))),obj(j),theID2(b(badID(i))),ID(j),theflux2(b(badID(i)))/flux(j)
                 if nfluxratio lt (1.0+delta2*(1/snr2(b(badID(i)))))^(-1) OR nfluxratio gt (1.0-delta2*(1/snr2(b(badID(i)))))^(-1) then begin
                    print,'Found extra bad flux.'
                    nbadflux(badID(i))=1
                    nbadlabel(badID(i))=ID(j)
                    badfluxes(badID(i),0)=theflux2(b(badID(i)))
                    badfluxes(badID(i),1)=flux(j)
                    abadflux(a(badID(i)))=1
                 endif
              endif
           endfor
        endfor
        endif

        print,nbadflux
        badID=where(nbadflux eq 1)
	badflux=where((zratio lt 1.0+delta AND zratio gt 1.0-delta) AND (theID(a) eq theID2(b)) AND (fluxratio lt (1.0+delta2*(1/snr(a)))^(-1) OR fluxratio gt (1.0-delta2*(1/snr(a)))^(-1)))
	badz=where(zratio gt 1.0+delta OR zratio lt 1.0-delta)        
	match2,object,object2,aa,bb
	match2,object,obj2,test_aa,test_bb
	match2,object2,obj,test_cc,test_dd	
	list1=n_elements(aa)
	list2=n_elements(bb)

;Write out table of failed matches.	
	openw,1,prefix+'matchextraction.dat'
        openw,2,prefix+'goodmatches.dat'
	; First objects which do not match
	for i=0,list1-1 do begin
		flag=zero
                cont_flag='-'
                dacomment1=' '
                dacomment2=' '
                OVER=0
                for q=0,lta3-1 do begin
                    if object(i) eq obj3a(q) and OVER eq 0 then begin
                            if comment3a(q) ne ' ' then begin
                               dacomment1=' 1:'+comment3a(q)
                               print,object(i),dacomment1
                               OVER=1
                            endif
                    endif
                endfor
		if aa(i) eq -1 then begin
                   if ltt4 gt 0 then begin
                      OVER=0
                      for q=0,ltt4-1 do begin
                         if object(i) eq obj4(q) and OVER eq 0 then begin
                            flag=three
                            if comment4(q) ne ' ' then begin
                               dacomment2=' 2:'+comment4(q)
                               OVER=1
                            endif
                         endif
                      endfor
                      print,object(i),flag,dacomment1+dacomment2
                   endif
                   if test_aa(i) gt -1 then begin
                      if theID(i) eq ID2(test_aa(i)) AND flux2(test_aa(i))/fluxE2(test_aa(i)) le snr_limit then flag=two
			;	print, 'Hallo',obj2(test_aa(i)),object(i),snr(i),flux2(test_aa(i))/fluxE2(test_aa(i))	
			
                   endif
                   if thecont(i) eq 2 then cont_flag='c'
                   if thecont(i) eq 2 AND flag eq 3 then cont_flag='C'
	           printf,1,object(i),one,flag,zero,cont_flag,thez(i),na,theflux(i),na,theID(i),dacomment1+dacomment2,FORMAT='(I6,3I6,A6,2F16.6,2E16.6,A16,A-100)'
                endif else begin
                        zratio2=thez(i)/thez2(aa(i))
                        fluxratio2=theflux(i)/theflux2(aa(i))
                        if (theID(i) ne theID2(aa(i))) then flag=three
                        if (zratio2 lt 1.0+delta AND zratio2 gt 1.0-delta) AND (((theID(i) eq theID2(aa(i))) AND $
                           (fluxratio2 gt (1.0+delta2*(1/snr(i)))^(-1) AND fluxratio2 lt (1.0-delta2*(1/snr(i)))^(-1))) OR $
                           ((theID(i) ne theID2(aa(i))) AND abadflux(i) ne 1)) then begin
                            print,aa(i),i,object(i),object2(aa(i)),thez(i),thez2(aa(i)),theflux(i),theflux2(aa(i))
                            if ((theID(i) ne theID2(aa(i))) AND abadflux(i) ne 1) then print,'WARNING Will! What is this crap?'
                            if thecont(i) eq 2 OR thecont2(aa(i)) eq 2 then cont_flag='c'
                            if thecont(i) eq 2 AND thecont2(aa(i)) eq 2 then cont_flag='C'
                     
                        printf,2,object(i),one,one,flag,cont_flag,thez(i),thez2(aa(i)),theflux(i),theflux2(aa(i)),theID(i),dacomment1+dacomment2,FORMAT='(I6,3I6,A6,2F16.6,2E16.6,A16,A-100)'
                        endif
                endelse
	endfor
	for i=0,list2-1 do begin
		flag=zero
                cont_flag='-'
                dacomment1=' '
                dacomment2=' '
                OVER=0
                for q=0,lta4-1 do begin
                    if object2(i) eq obj4a(q) and OVER eq 0 then begin
                            if comment4a(q) ne ' ' then begin
                               dacomment2=' 2:'+comment4a(q)
                               print,object2(i),dacomment2
                               OVER=1
                            endif
                    endif
                endfor
                if bb(i) eq -1 then begin
                   if ltt3 gt 0 then begin
                      OVER=0
                      for q=0,ltt3-1 do begin
                         if object2(i) eq obj3(q) and OVER eq 0 then begin
                            flag=three
                            if comment3(q) ne ' ' then begin
                               dacomment1=' 1:'+comment3(q)
                               OVER=1
                            endif
                         endif
                      endfor
                   endif
			if test_cc(i) gt -1 then begin
                             if theID2(i) eq ID(test_cc(i)) AND flux(test_cc(i))/fluxE(test_cc(i)) le snr_limit then flag=two
	;			print, 'Hallo',obj(test_cc(i)),object2(i),snr2(i),flux(test_cc(i))/fluxE(test_cc(i))
                  
                        endif
                        if thecont2(i) eq 2 then cont_flag='c'
                        if thecont2(i) eq 2 AND flag eq 3 then cont_flag='C'
			printf,1,object2(i),flag,one,zero,cont_flag,na,thez2(i),na,theflux2(i),theID2(i),dacomment1+dacomment2,FORMAT='(I6,3I6,A6,2F16.6,2E16.6,A16,A-100)'
                endif
        endfor		

	;Now objects with mismatched redshifts
	badlth=n_elements(badz)
        if badz(0) ne -1 then begin
             ;print,badz
             for i=0,badlth-1 do begin
                dacomment1=' '
                dacomment2=' '
                OVER=0
                for q=0,lta3-1 do begin
                    if object(a(badz(i))) eq obj3a(q) and OVER eq 0 then begin
                            if comment3a(q) ne ' ' then begin
                               dacomment1=' 1:'+comment3a(q)
                               print,object(a(badz(i))),dacomment1
                               OVER=1
                            endif
                    endif
                endfor
                OVER=0
                for q=0,lta4-1 do begin
                    if object2(b(badz(i))) eq obj4a(q) and OVER eq 0 then begin
                            if comment4a(q) ne ' ' then begin
                               dacomment2=' 2:'+comment4a(q)
                               print,object2(b(badz(i))),dacomment2
                               OVER=1
                            endif
                    endif
                endfor

                cont_flag='_'
                if thecont(a(badz(i))) eq 2 OR thecont2(b(badz(i))) eq 2 then cont_flag='c'
                if thecont(a(badz(i))) eq 2 AND thecont2(b(badz(i))) eq 2 then cont_flag='C'
		printf,1,object(a(badz(i))),one,one,one,cont_flag,thez(a(badz(i))),thez2(b(badz(i))),theflux(a(badz(i))),theflux2(b(badz(i))),theID(a(badz(i))),theID2(b(badz(i))),$
                       dacomment1+dacomment2,FORMAT='(I6,3I6,A6,2F16.6,2E16.6,2A16,A-100)'
             endfor
        endif

	;Now objects with mismatched fluxes
	badlth2=n_elements(badflux)
        if badflux(0) ne -1 then begin
             ;print,badflux
             for i=0,badlth2-1 do begin
                dacomment1=' '
                dacomment2=' '
                OVER=0
                for q=0,lta3-1 do begin
                    if object(a(badflux(i))) eq obj3a(q) and OVER eq 0 then begin
                            if comment3a(q) ne ' ' then begin
                               dacomment1=' 1:'+comment3a(q)
                               print,object(a(badflux(i))),dacomment1
                               OVER=1
                            endif
                    endif
                endfor
                OVER=0
                for q=0,lta4-1 do begin
                    if object2(b(badflux(i))) eq obj4a(q) and OVER eq 0 then begin
                            if comment4a(q) ne ' ' then begin
                               dacomment2=' 2:'+comment4a(q)
                               print,object2(b(badflux(i))),dacomment2
                               OVER=1
                            endif
                    endif
                endfor
                cont_flag='_'
                if thecont(a(badflux(i))) eq 2 OR thecont2(b(badflux(i))) eq 2 then cont_flag='c'
                if thecont(a(badflux(i))) eq 2 AND thecont2(b(badflux(i))) eq 2 then cont_flag='C'
                printf,1,object(a(badflux(i))),one,one,two,cont_flag,thez(a(badflux(i))),thez2(b(badflux(i))),theflux(a(badflux(i))),theflux2(b(badflux(i))),theID(a(badflux(i))),$
                       dacomment1+dacomment2,FORMAT='(I6,3I6,A6,2F16.6,2E16.6,A16,A-100)'
             endfor
          endif
                                ; Now objects with matching redshifts,
                                ; but mismatched line IDs, i.e. a
                                ; bright line has been missed (or mis-measured).
                                ; Checked for flux
                                ; mismatches. Included if there is one.
                              
	badlth3=n_elements(badID)
        if badID(0) ne -1 then begin
             ;print,badID,'badID'
             for i=0,badlth3-1 do begin
                printf,1,object(a(badID(i))),one,one,three,thez(a(badID(i))),thez2(b(badID(i))),badfluxes(badID(i),0),badfluxes(badID(i),1),nbadlabel(badID(i)),$
                       FORMAT='(I6,3I6,2F16.6,2E16.6,A16)'
             endfor
        endif
	close,1
        close,2
	readcol,prefix+'matchextraction.dat',number,user1,user2,flag,format='i,i,i,i'
	readcol,file5,fd,grism,number2,wave,numpix,sigma,flag2,format='i,a,i,f,i,f,i'
	match2,number2,number,aa,bb
	lth3=n_elements(aa)
	openw,1,prefix+'lines2.dat'
	for i=0,lth3-1 do begin
		if aa(i) ne -1 then begin
			printf,1,fd(i),grism(i),number2(i),wave(i),numpix(i),sigma(i),flag2(i),format='(I3,A9,I6,F18.6,I6,F14.6,I10)'
		endif
	endfor
	close,1

	print,'Objects in List1 with SNR > 3:',n_elements(object)
	print,'Objects in List2 with SNR > 3:',n_elements(object2)
	print,'Number of Multi-line Objects, List1:',n_elements(where(match eq 1))
	print,'Number of Multi-line Objects, List2:',n_elements(where(match2 eq 1))
end 	
