from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits as f
from astropy import units as u
from specutils.spectra import Spectrum1D
from astropy.modeling.models import Polynomial1D
from specutils.fitting.continuum import fit_continuum
from scipy.optimize import curve_fit

d346 = f.getdata('D346.fits')
d348 = f.getdata('D348.fits')
d350 = f.getdata('D350.fits')

header_d346 = f.getheader('D346.fits')
header_d348 = f.getheader('D348.fits')
header_d350 = f.getheader('D350.fits')

def mygaussian(x, amp, center, width):
    return 1.0 - amp*np.exp(-(x-center)**2/(2.0*width**2))

def get_desried_wavelength(file, header, flag=False):
    pixel_size = np.arange(file.shape[0])
    wl = header['CRVAL1'] + pixel_size*header['CD1_1']
    idx, = np.where(wl <= 6850 )
    wl_range = wl[idx]
    file = file[idx]
    spectrum = Spectrum1D(flux=file*u.adu, spectral_axis = wl_range*u.AA )
    p = Polynomial1D(degree=7)
    fit_res = fit_continuum(spectrum, model=p)
    fitted_spec = fit_res(wl_range * u.AA)
    
    plt.plot(spectrum.spectral_axis, spectrum.flux, 'r-')
    plt.grid(True)
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (ADU)')
    plt.plot(wl_range, fitted_spec, 'b-', label='D346')
    plt.legend(('D346 Spectrum', 'Fit',))
    
    # Let's divide the spectrum by the continuum model
    norm_flux = spectrum.flux / fitted_spec
    
    init_params = [0.3, 6562.7, 1.0]
    if flag:
        init_params = [0.3, 6717.6, 1.0]
    popt, pcov = curve_fit(mygaussian, wl_range, norm_flux, init_params)
    best_fit = mygaussian(wl_range, *popt)
    
    plt.figure()
    plt.grid(True)
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (ADU)')
    plt.plot(wl_range, norm_flux, 'r-', label='D346')
    plt.plot(wl_range, best_fit, 'b-', label='Continuum Fit')
    plt.legend(('D346', 'Continuum Fit',))
    plt.xlim(6530,6600)
    
    print(f'BFIT: line center = {popt[1]:.3f} A')
    print(f'BFIT: amplitude = {popt[0]:.3f}')
    print(f'BFIT: line width = {popt[2]:.3f} A')
    print(f'BFIT: equivalent width = {popt[0]*popt[2]:.3f} A')
    
    return norm_flux, wl_range, file

def measure_EW_area(flux, wavelength, wvl_lower, wvl_upper):
    ''' This function measures the equivalent width of a spectral line.
    
    Inputs:
    flux, wavelength -- spectrum to be measured
    wvl_lower, wvl_upper -- wavelength limits for the spectral line
    
    Outputs:
    EW -- equivalent width
    '''
    
    # Crop the spectrum to be within the wavelength limits
    idx = np.where((wavelength > wvl_lower) & (wavelength < wvl_upper))[0]
    print('idx', idx)
    
    # Define the wavelength difference between two adjacent points in the spectrum
    flux2 = [0.56, 1]
    wavelength2 = [wavelength[114], wavelength[116]]    
    EW = np.trapz(flux2, wavelength2)
    return EW

norm_flux_d346, wl_d346, d346 = get_desried_wavelength(d346, header_d346)
norm_flux_d348, wl_d348, d348 = get_desried_wavelength(d348, header_d348)
norm_flux_d350, wl_d350, d350 = get_desried_wavelength(d350, header_d350)
norm_flux_d346, wl_d346, d346 = get_desried_wavelength(d346, header_d346, flag=True)
norm_flux_d348, wl_d348, d348 = get_desried_wavelength(d348, header_d348, flag=True)
norm_flux_d350, wl_d350, d350 = get_desried_wavelength(d350, header_d350, flag=True)
