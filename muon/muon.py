import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2

def curve_fit_custom(f, tdata, ydata, p0=None, sigma=None, **kw):
    """
    Pass all arguments to curve_fit, which uses non-linear least squares
    to fit a function, f, to data.  Calculate the uncertaities in the
    fit parameters from the covariance matrix.
    """
    popt, pcov = curve_fit(f, tdata, ydata, p0, sigma, **kw)

    if sigma is None:
        chi2 = sum(((f(tdata,*popt)-ydata))**2)
    else:
        chi2 = sum(((f(tdata,*popt)-ydata)/sigma)**2)
    dof = len(ydata) - len(popt)
    rchi2 = chi2/dof
    print ('results of general_fit:')
    print ('   chi squared = ', chi2)
    print ('   degrees of freedom = ', dof)
    print ('   reduced chi squared = ', rchi2)
    
    return popt, chi2, dof

FILE_NAME1 = "./muon_148mV_corrected.txt"
FILE_NAME2 = "./muon_190mV_corrected.txt"
FILE_NAME3 = "./muon_260mV_corrected.txt"
FILE_NAME4 = "./muon_550mV_corrected.txt"

def model_func(x, tau, bg, amp):
    return bg + amp * np.exp(-1.0 * x / tau)

def slope(x, m, b):
    return m*x + b

def muon(filename, x, y, ratio, title):
    data = np.loadtxt(filename, usecols = (0))
    counts, tdata = np.histogram(data, bins='auto', density= True)
    tdata = (tdata[:-1] + tdata[1:])/2.
    sdev = np.sqrt(counts)/1000
    sdev[sdev == 0] = 9.31648059e-07
    mean = np.mean(data) 
    std = np.std(data)
    var = np.var(data)
    
    popt, chi, dof = curve_fit_custom(model_func, tdata, counts, p0 = [200, 10, .001], sigma=sdev)
    modelpoints = model_func(tdata, popt[0], popt[1], popt[2])
    
    plt.figure(figsize=(10,10))
    plt.grid()
    plt.plot(tdata, modelpoints, color = 'b', label = 'Exponential Fit Curve')
    plt.hist(data, bins = 'auto', density = True, color = "yellow", edgecolor='black', label = 'Histogram')
    plt.errorbar(tdata, counts, sdev, fmt='ro', ecolor=None, elinewidth=None,
         capsize=3, barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False, label = 'Error Bars')
    plt.title('{}mV Discriminator Treshold'.format(title), size = 18)
    plt.xlabel('Muon Decay Time (ns)', size = 14)
    plt.ylabel('Normalized Events/Bins', size = 14)
    plt.legend(loc='upper center', shadow=True, prop={'size': 14})
    plt.text(x, y, r"""$\chi^2$/ DOF ratio $\approx$ {}
    mean $\approx$ {} ns
    mean decay rate $\approx$ {} ns
    normalized amplitude $\approx$ {} ns""".format(ratio, int(mean), int(popt[0]), int(1/popt[2])),
             style='italic', fontsize=16, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
    plt.show()
    return tdata, modelpoints
    
def muonLog(filename, bg):
    data = np.loadtxt(filename, usecols = (0))
    counts, tdata = np.histogram(data, bins='auto')
    counts = counts - bg
    counts = np.log(counts[:35])
    counts[counts == -np.inf] = 0
    tdata = tdata[:-1]
    tdata = tdata[:35]
    sdev = np.sqrt(counts)/10
    sdev[sdev == 0] = 1

    popt, chi, dof = curve_fit_custom(slope, tdata, counts, p0 = [-200, 1], sigma=sdev)
    modelpoints = slope(tdata, popt[0], popt[1])
    
    plt.figure(figsize=(10,10))
    plt.grid()
    plt.plot(tdata, modelpoints, color = 'b')
    plt.errorbar(tdata, counts, sdev, fmt='ro', ecolor=None, elinewidth=None,
         capsize=3, barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False)
    plt.legend(('Curve fit', 'Error bars'), loc='upper center', shadow=True, prop={'size': 14})
    plt.xlabel('Muon Decay Time (ns)', size = 14)
    plt.ylabel('Events/Bins (Log Scale)', size = 14)
    plt.title('550mV Discriminator Treshold', size = 18)
    plt.text(6000, 5, r"""mean decay rate $\approx$ {} ns""".format(int(1/popt[0])),
             style='italic', fontsize=16, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
    plt.show()
    
muon(FILE_NAME1, 8000, .00015, 11.6, 148)
muon(FILE_NAME2, 8000, .00023, 1.9, 190)
muon(FILE_NAME3, 8000, .00025, 1.6, 260)
muon(FILE_NAME4, 8000, .0003, 1.1, 550)

muonLog(FILE_NAME4, 0)
#muonLog(FILE_NAME1, 70)
#muonLog(FILE_NAME2, 8)
#muonLog(FILE_NAME3, 12)
