from astropy.io import fits as f
import numpy as np
import matplotlib.pyplot as plt
from fixbads import *
from astropy.visualization import simple_norm

###############################################################################
## Step 1: Preprocessing D34-D45 data fits files
###############################################################################

d46_bias_i705 = f.getdata('D46.fits')
d47_bias_i705 = f.getdata('D47.fits')
d48_bias_i705 = f.getdata('D48.fits')
d49_bias_i705 = f.getdata('D49.fits')
d50_bias_i705 = f.getdata('D50.fits')
d51_bias_i705 = f.getdata('D51.fits')
d52_bias_i705 = f.getdata('D52.fits')

biases = np.stack((d46_bias_i705, d47_bias_i705, d48_bias_i705, d49_bias_i705, d50_bias_i705, d51_bias_i705, d52_bias_i705))
bias = np.median(biases, axis=0)

d01_dome_b639 = f.getdata('D01.fits')
d02_dome_b639 = f.getdata('D02.fits')
d03_dome_b639 = f.getdata('D03.fits')
d04_dome_b639 = f.getdata('D04.fits')
d05_dome_b639 = f.getdata('D05.fits')

d06_dome_v641 = f.getdata('D06.fits')
d07_dome_v641 = f.getdata('D07.fits')
d08_dome_v641 = f.getdata('D08.fits')
d09_dome_v641 = f.getdata('D09.fits')
d10_dome_v641 = f.getdata('D10.fits')

d11_dome_i705 = f.getdata('D11.fits')
d12_dome_i705 = f.getdata('D12.fits')
d13_dome_i705 = f.getdata('D13.fits')
d14_dome_i705 = f.getdata('D14.fits')
d15_dome_i705 = f.getdata('D15.fits')

d01_dome_b639_wo_bias = d01_dome_b639 - bias
d02_dome_b639_wo_bias = d02_dome_b639 - bias
d03_dome_b639_wo_bias = d03_dome_b639 - bias
d04_dome_b639_wo_bias = d04_dome_b639 - bias
d05_dome_b639_wo_bias = d05_dome_b639 - bias

d06_dome_v641_wo_bias = d06_dome_v641 - bias
d07_dome_v641_wo_bias = d07_dome_v641 - bias
d08_dome_v641_wo_bias = d08_dome_v641 - bias
d09_dome_v641_wo_bias = d09_dome_v641 - bias
d10_dome_v641_wo_bias = d10_dome_v641 - bias

d11_dome_i705_wo_bias = d11_dome_i705 - bias
d12_dome_i705_wo_bias = d12_dome_i705 - bias
d13_dome_i705_wo_bias = d13_dome_i705 - bias
d14_dome_i705_wo_bias = d14_dome_i705 - bias
d15_dome_i705_wo_bias = d15_dome_i705 - bias

domes_b639 = np.stack((d01_dome_b639_wo_bias, d02_dome_b639_wo_bias,
                  d03_dome_b639_wo_bias, d04_dome_b639_wo_bias,
                  d05_dome_b639_wo_bias))

domes_v641 = np.stack((d06_dome_v641_wo_bias,
                  d07_dome_v641_wo_bias, d08_dome_v641_wo_bias,
                  d09_dome_v641_wo_bias, d10_dome_v641_wo_bias))

domes_i705 = np.stack((d11_dome_i705_wo_bias, d12_dome_i705_wo_bias,
                  d13_dome_i705_wo_bias, d14_dome_i705_wo_bias,
                  d15_dome_i705_wo_bias))

dome_b639 = np.median(domes_b639, axis=0)
dome_v641 = np.median(domes_v641, axis=0)
dome_i705 = np.median(domes_i705, axis=0)

# We need to normalize our flats before we use them for calibration, so that we're not losing information about flux in our science images.
normalized_flat_b639 = dome_b639 / np.median(dome_b639)
normalized_flat_v641 = dome_v641 / np.median(dome_v641)
normalized_flat_i705 = dome_i705 / np.median(dome_i705)

d43_target_b639 = f.getdata('D43.fits')
d44_target_v641 = f.getdata('D44.fits')
d45_target_i705 = f.getdata('D45.fits')

# bias subtraction
d43_target_b639_wo_bias = d43_target_b639 - bias
d44_target_v641_wo_bias = d44_target_v641 - bias
d45_target_i705_wo_bias = d45_target_i705 - bias

# flat correction: notice that the flat field correction here divides instead of subtracts since unlike bias and dark frames that are “additive” terms,
# while flat effect is multiplicative.
d43_target_b639_wo_bias = d43_target_b639_wo_bias / normalized_flat_b639
d44_target_v641_wo_bias = d44_target_v641_wo_bias / normalized_flat_v641
d45_target_i705_wo_bias = d45_target_i705_wo_bias / normalized_flat_i705

badmask = f.getdata('badpixel_mask.fits')
d43_target_b639_final = fixbads(d43_target_b639_wo_bias, badmask, fixneg=False)
d44_target_v641_final = fixbads(d44_target_v641_wo_bias, badmask, fixneg=False)
d45_target_i705_final = fixbads(d45_target_i705_wo_bias, badmask, fixneg=False)

plt.figure(figsize=(10,10))
norm = simple_norm(d43_target_b639_final, 'sqrt', percent = 97)
plt.imshow(d43_target_b639_final, vmin=200, vmax=700, norm = norm)
f.writeto('preprocessedimages/D43_preprocessed.fits', d43_target_b639_final, header = f.getheader('D43.fits'), overwrite = True)
plt.figure(figsize=(10,10))
norm = simple_norm(d44_target_v641_final, 'sqrt', percent = 97)
plt.imshow(d44_target_v641_final, vmin=200, vmax=700, norm = norm)
f.writeto('preprocessedimages/D44_preprocessed.fits', d44_target_v641_final, header = f.getheader('D44.fits'), overwrite = True)
plt.figure(figsize=(10,10))
norm = simple_norm(d45_target_i705_final, 'sqrt', percent = 97)
plt.imshow(d45_target_i705_final, vmin=200, vmax=700, norm = norm)
f.writeto('preprocessedimages/D45_preprocessed.fits', d45_target_i705_final, header = f.getheader('D45.fits'), overwrite = True)

d34_pg_b639 = f.getdata('D34.fits')
d35_pg_v641 = f.getdata('D35.fits')
d36_pg_i705 = f.getdata('D36.fits')
d37_pg_b639 = f.getdata('D37.fits')
d38_pg_v641 = f.getdata('D38.fits')
d39_pg_i705 = f.getdata('D39.fits')
d40_pg_b639 = f.getdata('D40.fits')
d41_pg_v641 = f.getdata('D41.fits')
d42_pg_i705 = f.getdata('D42.fits')

# bias subtraction
d34_pg_b639_wo_bias = d34_pg_b639 - bias
d35_pg_v641_wo_bias = d35_pg_v641 - bias
d36_pg_i705_wo_bias = d36_pg_i705 - bias
d37_pg_b639_wo_bias = d37_pg_b639 - bias
d38_pg_v641_wo_bias = d38_pg_v641 - bias
d39_pg_i705_wo_bias = d39_pg_i705 - bias
d40_pg_b639_wo_bias = d40_pg_b639 - bias
d41_pg_v641_wo_bias = d41_pg_v641 - bias
d42_pg_i705_wo_bias = d42_pg_i705 - bias

# flat correction
d34_pg_b639_wo_bias = d34_pg_b639_wo_bias / normalized_flat_b639
d35_pg_v641_wo_bias = d35_pg_v641_wo_bias / normalized_flat_v641
d36_pg_i705_wo_bias = d36_pg_i705_wo_bias / normalized_flat_i705
d37_pg_b639_wo_bias = d37_pg_b639_wo_bias / normalized_flat_b639
d38_pg_v641_wo_bias = d38_pg_v641_wo_bias / normalized_flat_v641
d39_pg_i705_wo_bias = d39_pg_i705_wo_bias / normalized_flat_i705
d40_pg_b639_wo_bias = d40_pg_b639_wo_bias / normalized_flat_b639
d41_pg_v641_wo_bias = d41_pg_v641_wo_bias / normalized_flat_v641
d42_pg_i705_wo_bias = d42_pg_i705_wo_bias / normalized_flat_i705

d34_pg_b639_final = fixbads(d34_pg_b639_wo_bias, badmask, fixneg=False)
d35_pg_v641_final = fixbads(d35_pg_v641_wo_bias, badmask, fixneg=False)
d36_pg_i705_final = fixbads(d36_pg_i705_wo_bias, badmask, fixneg=False)
d37_pg_b639_final = fixbads(d37_pg_b639_wo_bias, badmask, fixneg=False)
d38_pg_v641_final = fixbads(d38_pg_v641_wo_bias, badmask, fixneg=False)
d39_pg_i705_final = fixbads(d39_pg_i705_wo_bias, badmask, fixneg=False)
d40_pg_b639_final = fixbads(d40_pg_b639_wo_bias, badmask, fixneg=False)
d41_pg_v641_final = fixbads(d41_pg_v641_wo_bias, badmask, fixneg=False)
d42_pg_i705_final = fixbads(d42_pg_i705_wo_bias, badmask, fixneg=False)

###############################################################################
### Step 2: Photometry
###############################################################################

from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from photutils import centroid_com,centroid_sources

# centroid_sources returns back more refined positions of the star based off the positions you pass in.
def performCentroidSources(image, x_list, y_list):      
    xs, ys = centroid_sources(image, x_list, y_list, box_size = 21)
    return xs, ys

def performCircularApertureAndAnnulus(image, positions, r, r_in, r_out, fitsFile = '', undectected_mag = ''):
    src_ap = CircularAperture(positions, r = r)
    sky_ap = CircularAnnulus(positions, r_in = r_in, r_out = r_out)
    apertures = [src_ap, sky_ap]
    phot_table = aperture_photometry(image, apertures)
    phot_table['xcenter'].info.format='%.1f'
    phot_table['ycenter'].info.format='%.1f'
    phot_table['aperture_sum_0'].info.format='%.1f'
    phot_table['aperture_sum_1'].info.format='%.1f'
    
    bkg_mean = phot_table['aperture_sum_1']/sky_ap.area
    bkg_sum = bkg_mean * src_ap.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    exptime = f.getheader(fitsFile)['exptime']
    fluxes = final_sum / exptime
    if undectected_mag:
        # resetting the aperture sum since we already determined the magnitude therefore there is no need to compute the aperture_sum_0 or aperture_sum_1
        phot_table['aperture_sum_0'][3] = 0
        phot_table['aperture_sum_1'][3] = 0
        fluxes[3] = 1
    mags = -2.5*np.log10(fluxes) + 25.0
    mags.info.format='%.3f'
    if undectected_mag:
        mags[3] = undectected_mag               
    phot_table.add_column(mags, name='mags')
    return mags

def undectectedStar(image):
    undetected_ap = CircularAperture([601, 490], r = 6.3)
    src_mask = undetected_ap.to_mask(method='center')
    src_data = src_mask.multiply(image)
    one_sigma = np.std(src_data)
    three_sigma = 3.0*one_sigma
    upperlimit_mag = -2.5*np.log10(three_sigma) + 25.0
    return upperlimit_mag

# Image positions for tagets D43-D45
x_list_target = (674.6, 597.8, 341.6, 601)
y_list_target = (349.7, 297.1, 510.5, 490)
        
xs_d43, ys_d43 = performCentroidSources(d43_target_b639_final, x_list_target, y_list_target)
d43_upperlimit_mag = undectectedStar(d43_target_b639_final)
xs_d44, ys_d44 = performCentroidSources(d44_target_v641_final, x_list_target, y_list_target)
d44_upperlimit_mag = undectectedStar(d44_target_v641_final)
xs_d45, ys_d45 = performCentroidSources(d45_target_i705_final, x_list_target, y_list_target)
d45_upperlimit_mag = undectectedStar(d45_target_i705_final)

positions_d43 = list(zip(xs_d43, ys_d43))
positions_d44 = list(zip(xs_d44, ys_d44))
positions_d45 = list(zip(xs_d45, ys_d45))

d43_unc_mag_b = performCircularApertureAndAnnulus(d43_target_b639_final, positions_d43, 7.2, 9.5, 12.5, fitsFile = 'D43.fits', undectected_mag = d43_upperlimit_mag)
d44_unc_mag_v = performCircularApertureAndAnnulus(d44_target_v641_final, positions_d44, 7.2, 9.5, 12.5, fitsFile = 'D44.fits', undectected_mag = d44_upperlimit_mag)
d45_unc_mag_i = performCircularApertureAndAnnulus(d45_target_i705_final, positions_d45, 7.2, 9.5, 12.5, fitsFile = 'D45.fits', undectected_mag = d45_upperlimit_mag)

x_list_pg1323_bvi = (96.7, 731.6, 791.3, 757.2) # x positions for primary, star A, star B, star C
y_list_pg1323_bvi = (513.4, 254.4, 123.5, 689.4) # y positions for primary, star A, star B, star C

xs_d34_b, ys_d34_b = performCentroidSources(d34_pg_b639_final, x_list_pg1323_bvi, y_list_pg1323_bvi)
positions_d34_b = list(zip(xs_d34_b, ys_d34_b))
d34_unc_mag_b = performCircularApertureAndAnnulus(d34_pg_b639_final, positions_d34_b, 5.8, 8.5, 11.5, fitsFile = 'D34.fits')

xs_d34_v, ys_d34_v = performCentroidSources(d35_pg_v641_final, x_list_pg1323_bvi, y_list_pg1323_bvi)
positions_d35_v = list(zip(xs_d34_v, ys_d34_v))
d35_unc_mag_v = performCircularApertureAndAnnulus(d35_pg_v641_final, positions_d35_v, 5.8, 8.5, 11.5, fitsFile = 'D35.fits')

xs_d34_i, ys_d34_i = performCentroidSources(d36_pg_i705_final, x_list_pg1323_bvi, y_list_pg1323_bvi)
positions_d36_i = list(zip(xs_d34_i, ys_d34_i))
d36_unc_mag_i = performCircularApertureAndAnnulus(d36_pg_i705_final, positions_d36_i, 5.8, 12, 15, fitsFile = 'D36.fits')

x_list_pg1525_bvi = (463.9, 576.5, 636.3, 761.2, 491.1)# x positions for primary, star A, star B, star C
y_list_pg1525_bvi = (267.2, 399.6, 349.7, 780.9, 240.2)# y positions for primary, star A, star B, star C

xs_pg1525_b, ys_pg1525_b = performCentroidSources(d37_pg_b639_final, x_list_pg1525_bvi, y_list_pg1525_bvi)
positions_pg1525_b = list(zip(xs_pg1525_b, ys_pg1525_b))
d37_unc_mag_b = performCircularApertureAndAnnulus(d37_pg_b639_final, positions_pg1525_b, 5.8, 9, 12, fitsFile = 'D37.fits')

xs_pg1525_v, ys_pg1525_v = performCentroidSources(d38_pg_v641_final, x_list_pg1525_bvi, y_list_pg1525_bvi)
positions_pg1525_v = list(zip(xs_pg1525_v, ys_pg1525_v))
d38_unc_mag_v = performCircularApertureAndAnnulus(d38_pg_v641_final, positions_pg1525_v, 5.8, 9, 12, fitsFile = 'D38.fits')

xs_pg1525_i, ys_pg1525_i = performCentroidSources(d39_pg_i705_final, x_list_pg1525_bvi, y_list_pg1525_bvi)
positions_pg1525_i = list(zip(xs_pg1525_i, ys_pg1525_i))
d39_unc_mag_i = performCircularApertureAndAnnulus(d39_pg_i705_final, positions_pg1525_i, 5.8, 8.5, 11.5, fitsFile = 'D39.fits')

x_list_pg1633_bvi = (98.5, 220.4, 675.1, 919.1)# x positions for primary, star A, star B, star C
y_list_pg1633_bvi = (717.4, 733.4, 354.1, 337.4)# y positions for primary, star A, star B, star C)

xs_pg1633_b, ys_pg1633_b = performCentroidSources(d37_pg_b639_final, x_list_pg1633_bvi, y_list_pg1633_bvi)
positions_pg1633_b = list(zip(xs_pg1633_b, ys_pg1633_b))
d40_unc_mag_b = performCircularApertureAndAnnulus(d40_pg_b639_final, positions_pg1633_b, 5.8, 11, 14, fitsFile = 'D40.fits')

xs_pg1633_v, ys_pg1633_v = performCentroidSources(d38_pg_v641_final, x_list_pg1633_bvi, y_list_pg1633_bvi)
positions_pg1633_v = list(zip(xs_pg1633_v, ys_pg1633_v))
d41_unc_mag_v = performCircularApertureAndAnnulus(d41_pg_v641_final, positions_pg1633_v, 5.8, 12, 15, fitsFile = 'D41.fits')

x_list_pg1633_bvi_2 = (102.5, 223.4, 677.1, 921.0)# x positions for primary, star A, star B, star C
y_list_pg1633_bvi_2 = (720.4, 735.4, 357.1, 339.7)# y positions for primary, star A, star B, star C)

xs_pg1633_i, ys_pg1633_i = performCentroidSources(d39_pg_i705_final, x_list_pg1633_bvi_2, y_list_pg1633_bvi_2)
positions_pg1633_i = list(zip(xs_pg1633_i, ys_pg1633_i))
d42_unc_mag_i = performCircularApertureAndAnnulus(d42_pg_i705_final, positions_pg1633_i, 5.8, 12, 15, fitsFile = 'D42.fits')

###############################################################################
## Step 3: Standardization to obtain the coefficients
###############################################################################

# extinction coefficient derived
kv0_b, kv0_v, kv0_i = 0.283, 0.133, 0.047
kv1_b, kv1_v, kv1_i = 0.049, 0.032, 0.000

# airmass for each standard star system
airmass_pg1323 = 1.078
airmass_pg1525 = 1.145
airmass_pg1633 = 1.610

# actual magnitudes for the pg1323 standard stars
pg1323_v = np.array([13.481, 13.591, 13.406, 14.003])
pg1323_bv = np.array([-0.140, 0.393, 0.761, 0.707])
pg1323_vi = np.array([-0.127, 0.506, 0.833, 0.759])

#extinction corrected magnitudes for pg1323 system
d34_b_cm = d34_unc_mag_b - (kv1_b + kv0_b*pg1323_bv)*airmass_pg1323
d35_v_cm = d35_unc_mag_v - (kv1_v + kv0_v*pg1323_bv)*airmass_pg1323
d36_i_cm = d36_unc_mag_i - (kv1_i + kv0_i*pg1323_vi)*airmass_pg1323

pg1525_v = np.array([15.053, 13.509, 16.403, 13.530, 16.301])
pg1525_bv = np.array([-0.198, 0.757, 0.730, 1.109, 0.564])
pg1525_vi = np.array([-0.168, 0.869, 0.808, 1.104, 0.757])

#extinction corrected magnitudes for pg1525 system
d37_b_cm = d37_unc_mag_b - (kv1_b + kv0_b*pg1525_bv)*airmass_pg1525
d38_v_cm = d38_unc_mag_v - (kv1_v + kv0_v*pg1525_bv)*airmass_pg1525
d39_i_cm = d39_unc_mag_i - (kv1_i + kv0_i*pg1525_vi)*airmass_pg1525

pg1633_v = np.array([14.397, 15.256, 12.969, 13.229])
pg1633_bv = np.array([-0.192, 0.873, 1.081, 1.134])
pg1633_vi = np.array([-0.212, 1.015, 1.090, 1.138])

#extinction corrected magnitudes for pg1633 system
d40_b_cm = d40_unc_mag_b - (kv1_b + kv0_b*pg1633_bv)*airmass_pg1633
d41_v_cm = d41_unc_mag_v - (kv1_v + kv0_v*pg1633_bv)*airmass_pg1633
d42_i_cm = d42_unc_mag_i - (kv1_i + kv0_i*pg1633_vi)*airmass_pg1633

pg1323_cm_b = (pg1323_v + pg1323_bv) - d34_b_cm
pg1323_cm_v = pg1323_v - d35_v_cm
pg1323_cm_i = (pg1323_v + pg1323_vi) - d36_i_cm

pg1525_cm_b = (pg1525_v + pg1525_bv) - d37_b_cm
pg1525_cm_v = pg1525_v - d38_v_cm
pg1525_cm_i = (pg1525_v + pg1525_vi) - d39_i_cm

pg1633_cm_b = (pg1633_v + pg1633_bv) - d40_b_cm
pg1633_cm_v = pg1633_v - d41_v_cm
pg1633_cm_i = (pg1633_v + pg1633_vi) - d42_i_cm

bv_list = np.concatenate([pg1323_bv, pg1525_bv, pg1633_bv])
vi_list = np.concatenate([pg1323_vi, pg1525_vi, pg1633_vi])

# Obtaining the standard coefficients
cm_b = np.concatenate([pg1323_cm_b, pg1525_cm_b, pg1633_cm_b])
ab1, ab0 = np.polyfit(bv_list, cm_b, 1)

cm_v = np.concatenate([pg1323_cm_v, pg1525_cm_v, pg1633_cm_v])
av1, av0 = np.polyfit(bv_list, cm_v, 1)

cm_i = np.concatenate([pg1323_cm_i, pg1525_cm_i, pg1633_cm_i])
ai1, ai0 = np.polyfit(vi_list, cm_i, 1)

std_pg1323_b = d34_b_cm + ab1*pg1323_bv + ab0
std_pg1323_v = d35_v_cm + av1*pg1323_bv + av0
std_pg1323_i = d36_i_cm + ai1*pg1323_vi + ai0

std_pg1525_b = d37_b_cm + ab1*pg1525_bv + ab0
std_pg1525_v = d38_v_cm + av1*pg1525_bv + av0
std_pg1525_i = d39_i_cm + ai1*pg1525_vi + ai0

std_pg1633_b = d40_b_cm + ab1*pg1633_bv + ab0
std_pg1633_v = d41_v_cm + av1*pg1633_bv + av0
std_pg1633_i = d42_i_cm + ai1*pg1633_vi + ai0

###############################################################################
### Step 4: Derive instrumental magnitude for the target stars
###############################################################################

airmass = 1.67

bv_target = (d43_unc_mag_b - d44_unc_mag_v + (ab0 - av0))/(1-(ab1-av1))
vi_target = (d44_unc_mag_v - d45_unc_mag_i + (av0 - ai0))/(1-(av1-ai1))

d43_corrected_magnitude = d43_unc_mag_b - (kv1_b + kv0_b*bv_target)*airmass
d44_corrected_magnitude = d44_unc_mag_v - (kv1_v + kv0_v*bv_target)*airmass
d45_corrected_magnitude = d45_unc_mag_i - (kv1_i + kv0_i*vi_target)*airmass

std_d43 = d43_corrected_magnitude + ab1*bv_target + ab0
std_d44 = d44_corrected_magnitude + av1*bv_target + av0
std_d45 = d45_corrected_magnitude + ai1*vi_target + ai0
