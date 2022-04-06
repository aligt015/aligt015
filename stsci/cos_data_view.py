from astropy.table import Table
from astropy.io import fits

from astroquery.mast import Observations

from pathlib import Path

import matplotlib.pyplot as plt

import numpy as np

#################################### Exercise 2 #############################################

# data directory we are going to place in our downloaded files.
data_dir = Path('./data/')
data_dir.mkdir(exist_ok=True)

# First download the required data
nuv_downloads = Observations.download_products(Observations.get_product_list(Observations.query_criteria(obs_id='lbbd01020')), 
        download_dir=str(data_dir), extension='fits', mrp_only=True, cache=False)

# Get header information and print.
nuv_x1d_filepath = Path('./data/mastDownload/HST/lbbd01020/lbbd01020_x1dsum.fits')
nuv_asn_filepath = Path('./data/mastDownload/HST/lbbd01020/lbbd01020_asn.fits')
hdu = fits.open(nuv_x1d_filepath)
nuv_x1d_header0 = hdu[0].header['PROCTIME']
nuv_x1d_header1 = hdu[1].header['NUMFLASH']
print(nuv_x1d_header0, nuv_x1d_header1)

nuv_x1d_data = Table.read(nuv_x1d_filepath)
nuv_asn_data = Table.read(nuv_asn_filepath)
print(nuv_x1d_data)
print(nuv_asn_data)

#################################### Exercise 2 #############################################

# Download the required FUVB data
nuv_downloads = Observations.download_products(Observations.get_product_list(Observations.query_criteria(obs_id='lcxv13050')), 
        download_dir=str(data_dir), extension='fits', mrp_only=True, cache=False)

fuv_x1d_filepath = Path('./data/mastDownload/HST/lcxv13050/lcxv13050_x1dsum.fits')
fuv_asn_filepath = Path('./data/mastDownload/HST/lcxv13050/lcxv13050_asn.fits')
fuv_x1d_data = Table.read(fuv_x1d_filepath)

wvln, flux, segment = fuv_x1d_data[1]["WAVELENGTH", "FLUX", "SEGMENT"]
# We are asked to normalize the flux
flux = flux/np.nanmax(flux) # Or you could do flux /= np.nanmax(flux). Both are the same.

fig1, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=100)
ax.plot(wvln, flux, linestyle="-", linewidth=0.25, c='black', label=segment)
ax.set_title("Fig. 2.1\nSimple COS Segment FUVB Spectrum", size=20)
ax.set_xlabel('Wavelength [$\AA$]', size=12)
ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size=12)
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
