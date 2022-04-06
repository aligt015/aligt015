from astropy.table import Table
from astropy.io import fits

from astroquery.mast import Observations

from pathlib import Path

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
