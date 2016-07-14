;   Procedure name: CREATE_COLUMN
;
;   Porpoise:
;       Reform UCX 72 level HCHO column from bitpunch output into HCHO columns using
;       Sigma dimensions and creating the shape factors used in reformulating the AMF
;       in OMHCHO data files.
;
;   Requirements: gamap2, get_pedges (geoschem edges from pressure surface function)
;   
;   Example: create_column, '200502'
;       Does the following to create a file called hcho_200502.he5
;       pull out the hcho column and sundry data from the bitpunch output to create
;       a HDF-EOS5 file containing the following structure:
;         GROUP "/" {
;         DATASET "GC_UCX_HCHOColumns" {
;         DATATYPE  H5T_COMPOUND {
;           H5T_ARRAY { [91][144] H5T_IEEE_F64LE } "VCHCHO";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "NHCHO";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "SHAPEZ";
;           H5T_ARRAY { [91][144] H5T_IEEE_F64LE } "VCAIR";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "NAIR";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "SHAPESIGMA";
;           H5T_ARRAY { [91] H5T_IEEE_F32LE } "LATITUDE";
;           H5T_ARRAY { [144] H5T_IEEE_F32LE } "LONGITUDE";
;           H5T_ARRAY { [73][91][144] H5T_IEEE_F64LE } "PEDGES";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "PMIDS";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F64LE } "SIGMA";
;           H5T_ARRAY { [72][91][144] H5T_IEEE_F32LE } "BOXHEIGHTS";
;         }
;         DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
;         }}
;       With VCs in molecules/m2, densities in molecules/m3, Sigma and pressure dimensions
;       as well as apriori shape factors from Palmer 01
;

; Function to return linear or geometric midpoints from an array of pedges(or whatever)
function midpoints, pedges, geometric=geometric
    N = n_elements(pedges)
    inds=indgen(N-1)
    diffs=pedges[inds] - pedges[inds+1]
    mids= pedges[inds] + diffs/2.0
    if keyword_set(geometric) then $
        mids= sqrt(pedges[inds] * pedges[inds+1])

    return, mids
end


;;;;;;;;;;;;;;;;;;;;;;;;;
;;; METHOD HERE
;;;;;;;;;;;;;;;;;;;;;;;;;
pro create_column, yyyymm, PLOTS=PLOTS
    ;yyyymm='200501'

    ptn='./trac_avg*' + yyyymm +'*'
    file=file_search(ptn)
    print, 'reading ', file

    ; GC OUTPUT DIMENSIONS: 144, 91, 72)
    ; Read hcho ( molec/molec )
    flag= CTM_GET_DATABLOCK( ppbv, 'IJ-AVG-$', Tracer=20, File=file )
    ; Read psurf ( hPa )
    flag= CTM_GET_DATABLOCK( psurf, 'PEDGE-$', Tracer=1,  File=file )
    ; read number density of air( molec/m3 )
    flag= CTM_GET_DATABLOCK( Nair, 'BXHGHT-$', Tracer=4,  File=file )
    ; read the boxheights ( m )
    flag= CTM_GET_DATABLOCK( boxH, 'BXHGHT-$', Tracer=1,  File=file )
    
    ; 72 levels from geos 5 ucx, use gamap2 function to grab grid
    g5      = ctm_grid(ctm_type('GEOS5',res=2))
    lons    =g5.xmid
    lats    =g5.ymid
    nlons   =g5.imx
    nlats   =g5.jmx
    ctmpedges = g5.pedge
    ctmpmids  = g5.pmid
    
    pedges  = dblarr(nlons,nlats,73)
    pmids   = dblarr(nlons,nlats,72)
    Nhcho   = dblarr(nlons,nlats,72)
    S_z     = dblarr(nlons,nlats,72)
    S_sig   = dblarr(nlons,nlats,72)
    sigma   = dblarr(nlons,nlats,72)
    
    
    ; timing
    t1=systime(1)
    
    ;use psurf to get pressure edges, then find geometric midpoints
    for x=0, nlons-1 do begin
        for y=0,nlats-1 do begin
            ; Get pressure edges and geometric pressure midpoints
            pedges[x,y,*] = get_pedges(psurf[x,y])
            pmids[x,y,*] = midpoints(pedges[x,y,*],/geometric)
        endfor
    endfor
    
    ; Air density in molecules/m3
    Nair = double(Nair) 
    
    ; Total column AIR = vertical sum of (Nair * height)
    ; molecs/m2      = density(molecs/m3) * height(m)
    Omega_A = total( Nair * boxH , 3 )
    
    ; Density of HCHO in molecs/m3
    Nhcho   = double(ppbv)*double(Nair)*1d-9 ; molecs / m3
    
    ; TOTAL COLUMN HCHO = vertical sum of ( Nhcho * height)
    Omega_H   = total(Nhcho * boxH, 3) 
    
    ; Using Sigma coordinate: (1 at surface to 0 at toa)
    ; sigma = (P - P_t) / ( P_s - P_t )
    P_toa = pedges[*,*,-1]
    P_ttb = psurf - P_toa ; pressure difference from surface to TOA
    for zz=0,71 do begin
        Sigma[*,*,zz] = (pmids[*,*,zz] - P_toa) / P_ttb
    endfor
    
    ; normalized shape factors ( by definition from Palmer 01 )
    Omega_ratio = Omega_A / Omega_H 
    mixingratio= ppbv * 1d-9
    for ii=0, 71 do begin
        S_z[*,*,ii]  = Nhcho[*,*,ii] / Omega_H
        S_sig[*,*,ii] = Omega_ratio * mixingratio[*,*,ii]
    endfor
    
    ; VCHCHO: molecs/m2, number densities: molecs/m3
    structdata={VCHCHO:Omega_H, NHCHO:Nhcho, ShapeZ:S_z, $
                VCAir:Omega_A, NAir:Nair, ShapeSigma:S_sig, $
                latitude:lats, longitude:lons, PEdges:pedges, PMids:pmids, $
                Sigma:Sigma, boxheights:boxH}

    ; write a structure to hdf5 as a compound datatype
    fout = '/short/m19/jwg574/rundirs/geos5_2x25_UCX/trac_avg/hcho_'+yyyymm+'.he5'
    fid = H5F_CREATE(fout)
    
    datatype_id = H5T_IDL_CREATE(structdata)
    dataspace_id=H5S_CREATE_SIMPLE(1) ; not so simple..
    
    dataset_id = H5D_CREATE(fid,'GC_UCX_HCHOColumns',datatype_id,dataspace_id)
    H5D_WRITE, dataset_id, structdata
    
    H5S_CLOSE, dataspace_id
    H5T_CLOSE, datatype_id

    H5F_CLOSE, fid
    
    ; check out the file
    ;s=h5_parse('filename.h5',/read_data)
    ;help, /structure, s.datastr._data ; or s.GC_UCX_HCHOColumns._data?


    print, systime(1)-t1, 'Seconds to recalculate and save UCX HCHO profile'

    if keyword_set(PLOTS) then begin
        
        ;plot pedge exampels
        !p.multi = 0
        cgps_open, 'column_vertical_dimension_plot.png'
        cgdisplay, 1600,1400
        !x.margin=[10,10] ; this seems to do nothing
        ;!x.omargin=[5,15] ; this did nothing
        x=30
        y=30
        xarr=fltarr(73)+1.0
        ; First plot, direct comparison of pressure vs sigma levels
        ; Black = pressure mids, Blue = CTM pressure mids, Red = SIGMA
        cgplot, xarr, ctmpmids, charsize=2, /YLOG, title='Model Levels', $
            xtitle="", ytitle='hPa', yrange=[1000,0.1], thick=2, color='magenta',$
            psym=1, xrange=[0,3]
        cgplot, xarr, pmids[x,y,*] , color='black', thick=2, psym=1, /overplot
        ; Sigma levels on new axis (0 to 1)
        xarr=xarr+1.0
        ; set up axis
        cgAxis, YAxis=1, YRange=[1.,1e-4], ytitle="unitless", /Save
        cgplot, xarr, Sigma[x,y,*], thick=1, psym=1, color='red', /overplot
    
        cglegend, colors=['black','magenta','red'], $
                titles=['GC_pmids(hPa)','CTM_pmids(hPa)','Sigma(unitless)'],$
                charsize=2.0, alignment=0, length=0.0, thick=1.5,$
                location=[0.585,0.86], vspace=3.0, psym=[1,1,1]
    
        ; save out the plot
        cgps_close, /png
        
        ; Plot the ppb column on both dimensions
        cgps_open, 'column_ppb_plot.png'
        cgdisplay, 1600,1200
        !x.margin=[10,10]
        !x.omargin=[5,5]
        ; ppb on pressure and sigma levels
        cgplot, ppbv[x,y,*], pmids[x,y,*], charsize=2, /YLOG, $
            xtitle="(molecs HCHO)/(1e9 molecs air)", thick=3, color='black', $
            yrange=[1000,0.1]
        ; overplot dots onto data points
        cgplot, ppbv[x,y,*], pmids[x,y,*], thick=2, psym=1, color='black', /overplot
        
        cgAxis, YAxis=1, YRange=[1.,1e-4], YLog=1, ytitle="unitless", /Save
        
        cgplot, ppbv[x,y,*], Sigma[x,y,*], color='magenta', /overplot
        ; overplot dots onto data points
        cgplot, ppbv[x,y,*], Sigma[x,y,*], thick=2, psym=1, color='magenta', /overplot
        
        cglegend, colors=['black','magenta'], titles=['hPa','Sigma'], $
            charsize=2.0, alignment=0, length=0.05, thick=1.5, location=[0.5,0.86], vspace=2.0
        
        ; save out the plot:
        cgps_close,/png 
        
        ; plot density on both dimensions
        cgps_open, 'column_density_plot.png'
        cgdisplay, 1600,1200
        !x.margin=[10,10]
        !x.omargin=[5,5]
        ; density on pressure and sigma levels
        cgplot, Nhcho[x,y,*], pmids[x,y,*], charsize=2, /YLOG, /XLOG, $
            xtitle="(molecs HCHO)/m3", thick=3, color='black', $
            xrange=[1e10, 1e16] , yrange=[1000,0.1]
        ; overplot dots onto data points
        cgplot, Nhcho[x,y,*], pmids[x,y,*], thick=2, psym=1, color='black', /overplot
        
        cgAxis, YAxis=1, YRange=[1.,1e-4], YLog=1, ytitle="sigma", /Save
        
        cgplot, Nhcho[x,y,*], Sigma[x,y,*], color='magenta', /overplot
        ; overplot dots onto data points
        cgplot, Nhcho[x,y,*], Sigma[x,y,*], thick=2, psym=1, color='magenta', /overplot
        
        cglegend, colors=['black','magenta'], titles=['hPa','Sigma'], $
            charsize=2.0, alignment=0, length=0.05, thick=1.5, location=[0.5,0.86], vspace=2.0
        
        ; save out the plot:
        cgps_close,/png 
        
        ; plot some shape factors
        cgps_open, 'shape_factor_plot.png'
        cgdisplay, 1600,1200
        !x.margin=[10,10]
        !x.omargin=[5,5]
        ; density on pressure and sigma levels
        cgplot, S_sig[x,y,*], pmids[x,y,*], charsize=2, /YLOG, $
            xtitle="S_s", thick=3, color='black', $
            yrange=[1000,0.1] ;,xrange=[1e10, 1e16]
        ; overplot dots onto data points
        cgplot, S_sig[x,y,*], pmids[x,y,*], thick=2, psym=1, color='black', /overplot
        
        cgAxis, YAxis=1, YRange=[1.,1e-4], YLog=1, ytitle="sigma", /Save
        
        cgplot, S_sig[x,y,*], Sigma[x,y,*], color='magenta', /overplot
        ; overplot dots onto data points
        cgplot, S_sig[x,y,*], Sigma[x,y,*], thick=2, psym=1, color='magenta', /overplot
        
        cglegend, colors=['black','magenta'], titles=['hPa','Sigma'], $
            charsize=2.0, alignment=0, length=0.05, thick=1.5, location=[0.5,0.86], vspace=2.0
        
        ; save out the plot:
        cgps_close,/png 

        
        !P.multi=0
        !x.margin=0
    endif


end
