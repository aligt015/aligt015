This project shows an example of basic astronomical image reduction, photometry and standardization techniques. After preprocessing the images, the B, V and I magnitudes of the four science targets are derived in the standard system. The science target images are the D43, D44 and D45 fits files. The extinction coefficient and standard transformation coefficient are also obtained.

Also, in the code there is logic to simulate an undetected star. Keep in mind, when simulating an undetected star, an upper limit and standard deviation need to be calculated. First by estimating the standard deviation and then estimating the upper limit from that. A good approach to derive the standard deviation is by performing multiple aperture-photometry measurements in different blank regions of the image where you know there are not any sources and using that to derive the standard deviation.

NOTE: The circular apertures and sky annuli have different area (i.e., different number of pixels), and we need to take it into account. Also, in the code to display the aperture and sky region, I used unnecessary large aperture sizes and sky annuli radii only for visualization purposes.

Ways to improve this code is to use more sophisticated libarbries available on photoutils and astropy. For instance, take a look here:

https://photutils.readthedocs.io/en/stable/detection.html

https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html

https://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html
