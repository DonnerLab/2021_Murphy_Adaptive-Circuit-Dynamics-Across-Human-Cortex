

def run_fooof():
    
    from fooof import FOOOF
    import scipy.io as sio
    import numpy as np
    from os.path import join
    
    fname = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/agg/basepow/GROUP_AV_basepow.mat'
    
    mat = sio.loadmat(fname)   # load matlab file containing power spectra for all clusters
    nclust = mat['ga_spectra'].shape[0]   # get number of clusters 
    nperm = mat['nperm'][0,0]   # get number of permutations
    goodfreqs = (mat['freqs']<49) | ((mat['freqs']>51) & (mat['freqs']<99)) | (mat['freqs']>101) # indices for exlcuding line noise freqs
    cfreq = mat['freqs'][0,goodfreqs[0]]
    
    width_lim = ([(mat['freqs'][0,2]-mat['freqs'][0,1])*1.2, 12.0])  # define limits on width of gaussian periodic components
    fm = FOOOF(peak_width_limits=width_lim, max_n_peaks=3, background_mode='fixed')  # create FOOOF structure
    
    fparams = np.empty((2,nclust))  # prespecify empty variables
    fspectra = np.empty((cfreq.size,nclust))
    frsq = np.empty((nclust))
    
    ifparams = np.empty((2,nclust,nperm))
    ifrsq = np.empty((nclust,nperm))
    
    for c in range(0,nclust):
        cspec = mat['ga_spectra'][c,goodfreqs[0],0]   # exlcude line noise freqs
        fm.fit(cfreq,cspec)
        
        fparams[:,c] = fm.background_params_
        fspectra[:,c] = fm.fooofed_spectrum_
        frsq[c] = fm.r_squared_
        
        for i in range(0,nperm):
            cspec = mat['ga_spectra'][c,goodfreqs[0],i+1]   # exlcude line noise freqs
            fm.fit(cfreq,cspec)
            
            ifparams[:,c,i] = fm.background_params_
            ifrsq[c,i] = fm.r_squared_
    
    
    fname = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/agg/basepow/GROUP_AV_fooof_lin_fits.mat'
    sio.savemat(fname, dict(fparams=fparams, fspectra=fspectra, frsq=frsq, cfreq=cfreq, ifparams=ifparams, ifrsq=ifrsq))
    
