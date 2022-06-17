from astropy.table import Table
from astropy.io import fits

from astroquery.mast import Observations

from pathlib import Path

import matplotlib.pyplot as plt

import numpy as np

#################################### Exercise 1 #############################################

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

#################################### Exercise 2.1 #############################################

nRows = 5  # How many segments we wish to split the spectrum into
wvln, flux, fluxErr, segment = fuv_x1d_data[0]["WAVELENGTH", "FLUX", "ERROR", "SEGMENT"]

# We want only wavelength range from 1635 - 1675.
rangewv = np.where(np.logical_and(wvln >= 1635, wvln <= 1675)) # Or you could do: rangewv = (wvln > 1635) & (wvln < 1675)
rangewvln, rangeflux, rangeerr = wvln[rangewv], flux[rangewv], fluxErr[rangewv]
minx, maxx = min(rangewvln), max(rangewvln)
rangex = maxx - minx

fig = plt.figure(figsize=(14, 20))

for i in range(nRows):
    min_ = minx + i*rangex/nRows
    max_ = minx + (i+1)*rangex/nRows
    ax = plt.subplot(nRows, 1, i+1)

    if i == 0:  # A way to set Title, xlabel, and ylabel that will work independent of number of rows
        ax.set_title(f"Fig. 2.4{segment[-1]} \nSegment {segment} Spectrum split into segments", size=30)
    if i == nRows - 1:
        ax.set_xlabel("Wavelength [$\AA$]", size=30)
    if i == int(nRows/2):
        ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size=30)

    # Create the plot itself
    ax.errorbar(rangewvln, rangeflux, rangeerr, c=plt.cm.rainbow((i+1)/nRows), alpha=0.8, marker='.', markerfacecolor='k', markersize=2, mew=0)

    ax.set_xlim(min_, max_)
plt.tight_layout()
plt.show()

#################################### Exercise 3.1 #############################################

plt.figure(figsize=(12, 6))
for i in range(2):
    wvln, flux, fluxErr, dataQual, segment = nuv_x1d_data[i]["WAVELENGTH"], nuv_x1d_data[i]["FLUX"],\
        nuv_x1d_data[i]["ERROR"], nuv_x1d_data[i]["DQ"], nuv_x1d_data[i]["SEGMENT"]

    plt.plot(wvln[dataQual == 0], flux[dataQual == 0], linewidth=1, alpha=0.8,
             label=f"{segment} Cleaned")
    plt.scatter(wvln[dataQual != 0], flux[dataQual != 0], s=[12, 4][i], c='r', alpha=1,
                label=f"{segment} Pix out-of-bounds")
plt.legend(fontsize=18)
plt.title("Exercise 3.1", size=25)
plt.xlabel('Wavelength [$\AA$]', size=20)
plt.ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size=15)
plt.tight_layout()