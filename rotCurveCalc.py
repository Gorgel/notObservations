from __future__ import division
from astropy import constants as con
import scipy as sp
import matplotlib.pyplot as plt
import astropy as ap
import astropy.units as u
import numpy as np
import pyfits as fits

#initiates values
c = con.c
HaTheory = (656.28 * 10**-9 * u.meter)
HaTheory

#loads fits file
Science = fits.open('MegaScience.fits')

#calculates wavelengt scale
step = 1.26401159695618
xMin = 5803.90927853889
xMax = 5803.90927853889+step*2040
Lambda = np.linspace(xMin, xMax, 2040)

#laods flux values into variable y
y = Science[0].data

#iterates over the spectrums and plot them in the
#same window.
for spectrum in y:
    plt.figure(1)
    plt.plot(Lambda[580:620],spectrum[580:620])
    plt.title('H-alpha', fontsize=19)
    plt.xlabel("$Angstrom$", fontsize=18)
    plt.ylabel("$Flux$", fontsize=18)
    plt.grid()

#create empty array to store peak value wavelengths
observedLambda = np.array([])
for spectrum in y:
    spectrum = spectrum[580:620]     #cuts spectrum into interval around H-alpha
    HaPeak = max(spectrum)     #selects peak of H-alpha line
    index = sp.where(spectrum == HaPeak)     #selects index for peak value
    LambdaHaPeak = Lambda[index]     #uses above index to select corresponding wavelength
    observedLambda = sp.append(observedLambda, LambdaHaPeak)     #store H-alpha wavelength in array


observedLambda = observedLambda* u.Angstrom     #give wavelengt unit of Angstrom


Vrad = c*(observedLambda.si-HaTheory)/HaTheory     #calculates radial velocities using doppler formula
galRadius = np.linspace(0.1,1,len(y))     #create temporary radius

Vrad = Vrad/1e3     #convert  radial velocity to km/s

#plot rotation curve
plt.figure(2)
plt.plot(galRadius, Vrad.value, 'o', color='red')
plt.title('Rotation curve', fontsize=19)
plt.xlabel("relative radius", fontsize=18)
plt.ylabel("$Velocity [$m/s$]$", fontsize=18)
plt.grid()

plt.show()

