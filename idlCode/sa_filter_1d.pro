;+
;
;  Name                 SA_FILTER_1D
;
;  Purpose              Filter a 1-D array
;
;  Inputs               in - 1D array. Must be evenly spaced in time or height and
;                             missing values interpolated
;
;                       low_cutoff - lower cutoff limit 
;  
;                       high_cutoff - higher cutoff.
;                       
;                       Note: These should be specified in # of points of input array. For example
;                             if array is a time-series with each time-step 2minute, and
;                             you want a 10min - 30min bandpass filter, set: LOW_CUTOFF=5, HIGH_CUTOFF=15
;
;  Output               out - Bandpass filtered array
;
;  Comments             Future improvements could add keywords to specify LOW, HIGH or BANDPASS
;                         filter for FILTER procedure call.
;
;  History              19/01/15 Created.  
;
;-

PRO SA_FILTER_1D,in,low_cutoff,high_cutoff,out

    IF NOT KEYWORD_SET(LOW_CUTOFF) OR NOT KEYWORD_SET(HIGH_CUTOFF) THEN STOP
      
    ; Zero pad
    n=N_ELEMENTS(in) & base2 = FIX(ALOG(n)/ALOG(2) + 0.4999)   ; power of 2 nearest to N
    in2 = [in-MEAN(in,/NaN),FLTARR(2L^(base2 + 1) - n)]
    
    ; Filter
    out=FILTER(in2,LOW_CUTOFF,HIGH_CUTOFF,/TIME,/BANDPASS) 
    out=out[0:n-1]
    out=FLOAT(out) ; COMPLEX to float
            
    RETURN
    

