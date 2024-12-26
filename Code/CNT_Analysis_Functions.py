# Needed Libraries
import json
import sys
import matplotlib.cm as cm
import numpy as np
from math import factorial    
from matplotlib.patches import Rectangle
from matplotlib import patheffects as pe
from matplotlib import pyplot as plt

# Functions
def AddColorBar(cb, axis, fig, orientation='vertical', w=0.02, fontsize=8, label=None, ylabelpad=8, xlabelpad=0, tickvals=None, ticklabels=None, cb_x=None,cb_y=None):
    """
    Adds color bar to an image in a subplot
    Input:    cb          - colorbar data defined in an imshow call (matplotlib image)
              axis        - subplot on which to attach the colorbar (matplotlib subplot)
              fig         - figure in which the subplot is drawn (matplotlib figure)
    Optional: orientation - place color bar to left or bottom of image (str either 'horizontal' or 'vertical') Default: 'vertical'
              w           - colorbar width (fraction of figure width) Default: 0.02
              fontsize    - size of font (float) Default: 8
              label       - label for colorbar (str) Default: None
              ylabelpad   - labelpad between colorbar label and colorbar ticks when in 'vertical' configuration (float) Default: 8
              xlabelpad   - labelpad between colorbar label and colorbar ticks when in 'horizontal' configuration (float) Default: 0
              tickvals    - custom tick values for the colorbar (list/array of floats) Default: Set by matplotlib
              ticklabels  - custom tick label for the colorbar (list/array of str) Default: Set by matplotlib
    Output: None (matplotlib figure modified)
    """   
    xy=axis.get_position()
    if orientation=='horizontal': cbax=fig.add_axes([xy.x0,xy.y0-1.5*w,xy.width,w])
    else: cbax=fig.add_axes([xy.x1+w/2.,xy.y0,w,xy.height])
    if orientation=='horizontal': plt.colorbar(cb,cax=cbax,orientation='horizontal')
    else: plt.colorbar(cb,cax=cbax)
    cbax.tick_params(labelsize=fontsize,length=2,pad=0.5)
    if orientation=='horizontal':
        if cb_x is None: cb_x=0.5
        if cb_y is None: cb_y=0.0
        if label is not None: cbax.set_xlabel(label,labelpad=xlabelpad,fontsize=fontsize,x=cb_x,y=cb_y)    
        if tickvals is not None: cbax.set_xticks(tickvals)
        if ticklabels is not None: cbax.set_xticklabels(ticklabels)
    else:
        if label is not None: cbax.set_ylabel(label,labelpad=ylabelpad,fontsize=fontsize,rotation=270)
        if tickvals is not None: cbax.set_yticks(tickvals)
        if ticklabels is not None: cbax.set_xticklabels(ticklabels)
    return
    
def AddScaleBar(axis,scale,cal,w=None,xy=None,units='nm',fontsize=8,text=True):
    """
    Adds scale bar to an image in a subplot
    Input:    axis     - subplot axis for added scalebar (matplotlib axis)
              scale    - desired length of scalebar in desired unit (int)
              cal      - calibration of image in desired units/pixel (float)
    Optional: w        - width of scale bar in pixels (float) Default: 1/50th of image height
              xy       - coordinates of top left corner of scale bar (tuple) NOTE: coordinates of bar not text
                          Default: None (puts bar w pixels in each dimension away from bottom left corner)
              units    - desired units (str) Default: 'nm'
              fontsize - size of font (float) Default: 8
              text     - Display actual text of scale bar (bool) Default: True
    Output: fig - displayed image (matplotlib figure)
    """   
    imdim=axis.get_images()[0].get_extent()
    dim=[int(imdim[2]-imdim[3]),int(imdim[1]-imdim[-0])]
    if not w: w=dim[0]/50.
    if not xy: xy=(w,dim[0]-2*w)
    axis.add_patch(Rectangle(xy,scale/cal,w,fc='w',ec='k',lw=0.5))
    if text: 
        txt=axis.text(xy[0]+0.5*scale/cal,xy[1],str(scale)+' '+units,fontweight='bold',color='w',fontsize=fontsize,ha='center',va='bottom')
        txt.set_path_effects([pe.withStroke(linewidth=1,foreground='k')])
    return

def Display_SI(dat,disp,cal=1.,scale=None,en=None,cmap=cm.gist_gray,color='k',dpi=200,figsize=(9,3),fs=8):
    """
    Displays the spectral average (false bright field) of a spectrum image
    
    Input:    dat     - spectrum image data (numpy array)
    Optional: cal     - calibration of image if filename loaded, if data input directly (float)
              scale   - scale bar (float or int)
                         Default: None (calculates 1/5 of image width)
              cmap    - color map used for image display (python colormap from matplotlib) Default: Greys_r
              color   - color for spectrum plots (python color from matplotlib) Default: 'k'
              dpi     - dots per inch (int) Default: 200
              figsize - figure size in inches (tup): Default (4,3)
    Output:   none
    """    

    dim=dat.shape
    im=np.average(dat,axis=2)
    spec=np.average(dat,axis=(0,1))
    if en is None: en=GetUnalignedEnergyAxis(spec,disp)
    fwhm=Get_FWFM(en,spec,0.5)
    
    fig = plt.figure(figsize=figsize,dpi=dpi)
    relative_widths=[dim[1]/dim[0],1,1]
    gs = plt.GridSpec(ncols=3, nrows=1, width_ratios=relative_widths)
    ax=np.array([fig.add_subplot(gs[i]) for i in range(3)])
    
    plt.setp(ax[0],xticks=[],yticks=[])
    ax[0].imshow(im,cmap=cmap)
    if not scale: scale=int(dim[0]*cal/5.) 
    AddScaleBar(ax[0],scale,cal,fontsize=fs)
    
    plt.setp(ax[1:],xlabel='Energy Loss (eV)',ylabel='EELS Intensity (counts)')
    ax[1].plot(en,spec,color=color)
    ax[1].set_xlim(-3*fwhm[2],3*fwhm[2])
    ax[1].text(0.99,0.99,'FWHM\n'+str(round(fwhm[2]*1000,1))+' meV', va='top', ha='right', transform=ax[1].transAxes)
    ax[1].axvline(color='k',ls='--')
    ax[1].set_title('ZLP')
    
    plt.setp(ax[2],yscale='log')
    ax[2].plot(en,spec,color=color)
    ax[2].set_title('Spectrum')
    return

def GetCalibratedSpectra(dat,ens,en_lo=None,en_hi=None):
    """
    Creates a new dataset where all spectra in an arbitrarily sized sequence of EEL spectra (can be single spectrum, recording, linescan, spectrum image and 2D spectrum image) are aligned to a single energy axis. NOTE: Unlike previous iterations of this code this does not interpolate, it cuts elements off the edge to get everything aligned.
    
    Input:    dat   - data (numpy array)
              ens   - calibrated energy axes (numpy array)
    Optional: en_lo - desired lowest energy of the output data (float)
                        Default: None (finds highest minimum energy available for all spectra)
              en_hi - desired highest energy of cropped data (float)
                        Default: None (finds lowest maximum energy available for all spectra)
    Output:   E     - calibrated energy axis (numpy array)
              Specs - calibrated spectrum image (numpy array)
    """   

    print("Aligning Spectra to Calibrated Energy Axes")
    dim=dat.shape
    dat_flat=dat.reshape(-1,dim[-1])
    ens_flat=ens.reshape(-1,dim[-1])
    m=np.amax(ens_flat[:,0]);M=np.amin(ens_flat[:,-1])
    disp=np.mean(np.diff(ens_flat,axis=1))
    if not en_lo or en_lo<m: en_lo=m
    if not en_hi or en_hi>M: en_hi=M
    N=int((en_hi-en_lo)/disp)-1
    E=np.linspace(en_lo+disp,en_hi-disp,N)
    Specs,counter,total=[],0,np.prod(dim[:-1])
    for en,spec in zip(ens_flat,dat_flat):
        Specs.append(spec[Get_i(en,en_lo):Get_i(en,en_lo)+N])
        counter+=1;WritePercentage(counter,total)
    Specs=np.array(Specs).reshape((*dim[:-1],E.shape[0]))
    Spec_Av=np.average(Specs,axis=tuple(range(Specs.ndim-1)))
    E=E-Get_FWFM(E,Spec_Av,0.5,fitpoly=True)[-1]  
    print(' (FINISHED!)\n')
    return E,Specs

def GetExponentialFit_2R(en, spec, n, fran_lo, fran_hi, visualize = True):
    """
    Fits an exponential to an EEL spectrum over two energy ranges
    
    Input:    en       - energy axis of EEL spectrum (numpy array)
              spec     - EELS axis of EEL spectrum (numpy array)
              n        - order of exponential polynomial used for fit (int)
              fran_lo  - lower energy fit range  (float)
              fran_hi  - upper energy fit range (float)
    Optional: viualize - Plot output (boolean): Default - True
    Output:   en_out   - energy axis of the fitted region (numpy array)
              spec_out - EELS intensity in fitted region (numpy array)    
    """   
    
    i11 = Get_i(en, fran_lo[0]); i12 = Get_i(en, fran_lo[1])
    i21 = Get_i(en, fran_hi[0]); i22 = Get_i(en, fran_hi[1])
    en_fit   = np.append(en[i11:i12],   en[i21:i22])
    spec_fit = np.append(spec[i11:i12], spec[i21:i22])
    en_fit   = en_fit[np.where(spec_fit > 0)]
    spec_fit = spec_fit[np.where(spec_fit > 0)]
    spec_fit_log = np.log(spec_fit)
    fparams = np.polyfit(en_fit, spec_fit_log, n)
    en_out = en[i11:i22]
    bg_f = np.poly1d(fparams)
    bg = np.exp(bg_f(en_out))   
    spec_out = spec[i11:i22] - bg
    if visualize:
        f,a1=plt.subplots(1,1,figsize=(3.5,2.5),dpi=200)
        plt.setp(a1,xlabel='Energy Loss (meV)',ylabel='vEELS Intensity (counts)')
        
        vis_ran = np.array([fran_lo[0], fran_hi[1]])
        i1_vis = Get_i(en, vis_ran[0] - np.ptp(vis_ran)/2.)
        i2_vis = Get_i(en, vis_ran[1] + np.ptp(vis_ran)/2.)
        en_vis = en[i1_vis:i2_vis]; spec_vis = spec[i1_vis:i2_vis]
        ym,yM = 0.9 * np.amin(spec_vis), 1.025 * np.amax(spec_vis)
        
        a1.plot(en_vis, spec_vis, color='k', label='Data')
        a1.fill_betweenx([ym,yM],*fran_lo,color='b',alpha=0.25,lw=0,label='Fit Region')
        a1.fill_betweenx([ym,yM],*fran_hi,color='b',alpha=0.25,lw=0)
        a1.plot(en_out, bg, color='b', ls='--', label='Bkg. Fit')
        a1.set_ylim(ym,yM)
        
        a2=a1.twinx()
        a2.plot(en_out, spec_out, color='b', label='Fitted Peak')
        a2.axhline(0,color='k',ls='--',lw=0.5)
        a2.set_ylabel('Bkg. Sub. vEELS Intensity (counts)',color='b',rotation=270,labelpad=7)
        a2.tick_params(axis='y', labelcolor='b', colors='b')
        a2.spines['right'].set_color('b')
        
        h1, l1 = a1.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
        a1.legend(h1+h2,l1+l2,handlelength=0.75,loc='upper right')   
    return en_out, spec_out

def GetFilteredSpectrum(spec, w, order, deriv=0, rate=1):
    """
    Savitzky Golay filters a spectrum for noise along a certain window
    
    Input:  spec   - spectrum to be filtered (numpy array)
            w      - size of window (int)
            order  - order of polynomial fit (int)
    Output: spec_f - filtered spectrum (numpy array)
    """    
    try:
        w = np.abs(int(w))
        order = np.abs(int(order))
    except ValueError:
        raise ValueError("w and order have to be of type int")
    if w % 2 != 1 or w < 1:
        raise TypeError("w size must be a positive odd number")
    if w < order + 2:
        raise TypeError("w is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (w -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = spec[0] - np.abs( spec[1:half_window+1][::-1] - spec[0] )
    lastvals = spec[-1] + np.abs(spec[-half_window-1:-1][::-1] - spec[-1])
    spec = np.concatenate((firstvals, spec, lastvals))
    spec_f = np.convolve( m[::-1], spec, mode='valid')
    return spec_f

def Get_FWFM(en, spec, frac, norm='max', specsig=None, fitpoly=False, poly_o=2, fit_w=2, fitvis=False):
    """
    Finds FWHM of ZLP and returns fractional-maxes and width in energy
    
    Input:    en      - energy loss axis (numpy array)
              spec    - EELS axis (numpy array)
              frac    - fractional max of ZLP to find (frac)
    Optional: norm    - 'max': normalizes by dividing by max (Default)
                        'full': normalizes such that max is 1 and min is 0
              specsig - sigma to gaussian blur the spectrum for fwhm measurement (float): Default: None
              fitpoly - fit a polynomial to peak tails to obtain subpixel fwhm value (boolean) Default: False
              poly_o  - order of polynomial to be fit (int) Default: 2
              fit_w   - width of fit region for fitpoly (int) Default: 2
    Output:   en_lo   - energy of lower half max in eV (float)
              en_hi   - energy of upper half max in eV (float)
              en_w    - FWHM of ZLP in eV (float)
    """    
    from scipy import ndimage,interpolate
    if norm=='full': s=NormArray(spec)
    else: s=spec/np.amax(spec)
    if specsig: s=ndimage.gaussian_filter1d(s,specsig)
    i_lo=Get_i(s[:np.argmax(s)],frac)
    i_hi=Get_i(s[np.argmax(s):],frac)+np.argmax(s)
    if fitpoly:
        w=int(fit_w/2.)
        fit_e_lo=en[i_lo-w:i_lo+w+1];fit_s_lo=s[i_lo-w:i_lo+w+1]
        fit_e_hi=en[i_hi-w:i_hi+w+1];fit_s_hi=s[i_hi-w:i_hi+w+1]
        fit_lo=np.polyfit(fit_e_lo,fit_s_lo,poly_o)
        fit_hi=np.polyfit(fit_e_hi,fit_s_hi,poly_o)
        roots_lo=np.roots([*fit_lo[:-1],fit_lo[-1]-frac])
        roots_hi=np.roots([*fit_hi[:-1],fit_hi[-1]-frac])
        en_lo=roots_lo[np.argmin(np.abs(roots_lo-np.mean(fit_e_lo)))]
        en_hi=roots_hi[np.argmin(np.abs(roots_hi-np.mean(fit_e_hi)))]
    else: en_lo=en[i_lo];en_hi=en[i_hi]
    en_w=en_hi-en_lo;en_c=(en_hi+en_lo)/2.
    if fitvis:
        f,a=plt.subplots(1,1,dpi=200)
        a.axhline(0,color='k',ls=':',lw=0.5)
        a.axhline(frac,color='k',ls='--',lw=0.5)
        a.axhline(1,color='k',lw=0.5)
        a.plot(en,s,color='k',lw=0,marker='x')
        if fitpoly:
            a.plot(fit_e_lo,np.polyval(fit_lo,fit_e_lo),color='r',lw=1)
            a.plot(fit_e_hi,np.polyval(fit_hi,fit_e_hi),color='b',lw=1)
        a.axvline(en_lo,color='r',ls='--')
        a.axvline(en_hi,color='b',ls='--')
        a.set_xlim(en_lo-2*en_w,en_hi+2*en_w)
    return en_lo,en_hi,en_w,en_c

def Get_i(arr,val):
    """
    Gives index of array closest to desired value
    
    Input:  arr - numpy array
            val - float
    Output: i   - index of closest value to desired float (int)
    """    
    
    i=np.argmin(np.abs(arr-val))
    return i

def GetSpectraEnergyAxes(dat, disp, smoothsig=None, spec_thresh=10, style='fwfm', fit_fwfm=False, frac_fwfm=0.5, subpix_w=30.,subpix_g0=None,subpix_bounds=None):
    """
    Creates a numpy array of calibrated energy loss axes for an arbitrarily sized sequence of EEL spectra 
    
    Input:    dat          - EELS data (numpy array)
              disp         - dispersion of acquisition in eV (float)
    Optional: smoothsig    - sigma for gaussian blur applied to ZLP to find center (float) Default: None
              spec_thresh  - threshold for total number of counts in spectrum before ignoring.
              style        - 'pixel': calibrates to pixel with max intensity of ZLP
                             'fwhm': calibrates to center of fwhm of ZLP (Default)
                             'subpixel': fits Lorentzian to the ZLP
              fit_fwhm     - Fit a polynomial to the HM area instead of pixel FWHM (bool) Default: False
              frac_fwfm    - Fractional maximum of the ZLP to fit (frac) Default: 0.5
              subpix_w     - Number of pixels around ZLP for subpixel fit (int) Default: 30
              subpix_p0    - Initial guess for Lorentzian parameters for subpixel fit (lorz parameters) Default: (spec max, spec argmax, subpix_w/2.)
              subpix_bounds - Bounds for subpixel fit (scipy.optimize.curve_fit style bounds) Default: None
    Output:   ens          - calibrated energy axis for each spectra in the EELS data (numpy array)
    """        
    
    dim=dat.shape
    dat_flat=dat.reshape(-1,dim[-1])
    dat_av=np.average(dat_flat,axis=0)
    en_init=np.arange(0,dat.shape[-1])*disp
    ens,counter,total,c=[],0,np.prod(dim[:-1]),en_init[np.argmax(dat_av)]
    print("Finding Calibrated Energy Axes")
    for spec in dat_flat:
        if np.amax(spec)<spec_thresh: c=0.; continue
        if smoothsig is not None: spec=ndimage.gaussian_filter1d(spec,smoothsig)
        if style=='pixel': c=en_init[np.argmax(spec)]
        if style=='fwfm': c=np.average(Get_FWFM(en_init,spec,frac_fwfm,fitpoly=fit_fwfm)[:2])
        if style=='subpixel': 
            w=int(subpix_w/2.);i1=np.argmax(spec)-w;i2=np.argmax(spec)+w+1
            if subpix_g0 is None: subpix_g0=np.amax(spec),en_init[np.argmax(spec)],w
            c=optimize.curve_fit(lorz,en_init[i1:i2],spec[i1:i2],*subpix_g0,bounds=subpix_bounds)[0][1]
        ens.append(en_init-c)
        counter+=1;WritePercentage(counter,total)
    ens=np.asarray(ens).reshape(dim)
    print('\n')
    return ens 

def GetUnalignedEnergyAxis(dat, disp, style='fwhm'):
    """
    Creates an energy axis based off of the dispersion and an unaligned arbitrarily sized sequence of EEL spectra (can be recording, linescan, spectrum image and 2D spectrum image) NOTE: If data is a single spectrum the code will just provide a pixel fit of the spectrum, use GetSpectraEnergyAxes instead.
    
    Input:    dat   - EELS data (numpy array)
              disp  - dispersion of acquisition in eV (float)
    Optional: style - 'pixel': calibrates to pixel with max intensity of ZLP
                       'fwhm': calibrates to center of fwhm of ZLP (Default)
    Output:   en    - calibrated energy axis for the average spectrum in the unaligned EELS data (numpy array)
    """  
    dim=dat.shape
    dat_flat=dat.reshape(-1,dim[-1])
    dat_av=np.average(dat_flat,axis=0)
    if style=='pixel': dat_zlp=np.argmax(dat_av)
    if style=='fwhm': dat_zlp=Get_FWFM(np.arange(dim[-1]),dat,0.5)[-1]
    en=(np.arange(0,dim[-1])-dat_zlp)*disp 
    return en

def LoadEELS(fname):
    """
    Loads an EELS Spectrum or Recording
    
    Input:  fname - path to either .npy or .json file that (string)
    Output: dat   - data (numpy array)
            disp  - energy dispersion in eV (float)
            mdat  - metadata (python dictionary)
    """    

    if fname.endswith('npy'): fname=fname[:-4]
    if fname.endswith('json'): fname=fname[:-5]
    dat=np.load(fname+'.npy')
    mdat= json.load(open(fname+'.json','r'))
    disp=mdat['spatial_calibrations'][-1]['scale']
    return dat,disp,mdat

def LoadEELS_SI(fname):
    """
    Loads an EEL Spectrum Image
    
    Input:  fname - path to either .npy or .json file that (string)
    Output: dat   - data (numpy array)
            cal   - pixel calibration in nm (float)
            disp  - energy dispersion in eV (float)
            mdat  - metadata (python dictionary)
    """    

    if fname.endswith('npy'): fname=fname[:-4]
    if fname.endswith('json'): fname=fname[:-5]
    dat=np.load(fname+'.npy')
    mdat= json.load(open(fname+'.json','r'))
    cal=mdat['spatial_calibrations'][0]['scale']
    disp=mdat['spatial_calibrations'][-1]['scale']
    return dat,cal,disp,mdat
    
def LoadImage(fname):
    """
    Loads a 2D STEM Image
    
    Input:  fname - path to either .npy or .json file that (string)
    Output: im    - image data (numpy array)
            cal   - pixel calibration in nm (float)
            mdat  - metadata (python dictionary)
    """    

    if fname.endswith('npy'): fname=fname[:-4]
    if fname.endswith('json'): fname=fname[:-5]
    dat=np.load(fname+'.npy')
    mdat= json.load(open(fname+'.json','r'))
    cal=mdat['spatial_calibrations'][0]['scale']
    return dat,cal,mdat

def NormArray(arr):
    """
    Converts an array of arbitrary values to a 0-1 range
    
    Input:  arr  - numpy array
    Output: narr - normalized numpy array
    """    
    
    arr = np.asarray(arr)
    M=np.amax(arr);m=np.amin(arr)
    narr = (arr-m)/(M-m)
    return narr

def WritePercentage(i,total):
    """
    Writes out percent of total based off of index
    
    Input:  i     - index in for loop (int)
            total - total indeces to run (int)
    Output: Doesn't output anything, some magic I don;t understand happens that causes it to overwrite the previous %
    """    
    if i<total:
        if (i%(total/100.))<1:
            sys.stdout.write('\r')   
            sys.stdout.write("%d%%" % int(np.round(100*(i/total),0)))
            sys.stdout.flush()
            return
    elif i==total:
        sys.stdout.write('\r')   
        sys.stdout.write('100%')
        sys.stdout.flush()        
        return
    else:
        sys.stdout.write('\r')   
        sys.stdout.write('Something went wrong')
        sys.stdout.flush()        
        return 
