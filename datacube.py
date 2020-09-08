# Data Cube movie - Brian Kent and Jeff Mangum, NRAO
# September 2015
import os, sys, string, matplotlib, aplpy

import matplotlib.pyplot as plt
import numpy as np
import pylab as py
from matplotlib import colors, cm
from astropy.io import fits

# File name, number of channels, and header information
# Data from: http://www.mpia.de/THINGS/Data_files/NGC_2841_NA_CUBE_THINGS.FITS
filename = 'NGC_2841_NA_CUBE_THINGS.FITS'
nchan = 128
hdulist = fits.open(filename)  # open a FITS file
prihdr = hdulist[0].header           # the primary HDU header

# Colormap I like to use
cdict_bkcmap={
    'red'  :  ((0., 0., 0.), (0.25,0.,0.), (0.5,1.,1.), (0.75,1.0,1.0),  (1., 1., 1.)),
    'green':  ((0., 0., 0.), (0.25,0.,0.), (0.5,0.,0.), (0.75,1.0,1.0),  (1., 1., 1.)),
    'blue' :  ((0., 0., 0.), (0.25,1.,1.), (0.5,0.,0.), (0.75,0.0,0.0),  (1., 1., 1.))
    }

bkcmap = matplotlib.colors.LinearSegmentedColormap('bkcmap', cdict_bkcmap,1024)

# FITS coordinates for marking the spectral box rectangle
raminbox = 462 - 11
ramaxbox = 462 + 11
decminbox = 575 - 11
decmaxbox = 575 + 11

# Create 1D arrays for an averaged spectrum and velocity array
img = fits.getdata(filename)
spectrum = [np.mean(img[0,i,decminbox:decmaxbox,raminbox:ramaxbox])*1000.0 for i in range(0,nchan)]
velocity = [(prihdr['CRVAL3'] + prihdr['CDELT3']*(i+prihdr['CRPIX3']))/1000.0 for i in range(0,nchan)]

#Create figure canvas
fig = plt.figure(facecolor='w', edgecolor='w', frameon=True, figsize=(6,7))

#Remove old PNG frames
print 'WARNING: Removing old *.png files from pngs directory...'
os.system('rm -rf pngs/*.png')

for i in range(0,nchan):
        
        #Clear plot
        plt.clf()
        
        #Load and customize the data cube slice
        ax1 = aplpy.FITSFigure(filename, dimensions=[0,1], slices=[i,0], figure=fig, subplot=[0.25,0.5,0.60,0.45])
        ax1.recenter(140.50924, 50.975394, width=0.125, height=0.125)  # degrees	
        ax1.show_colorscale(cmap=bkcmap, vmin=-.0001, vmax=0.00387265)
        ax1.add_grid()
        ax1.axis_labels.set_xtext(u'$\\mathrm{Right}$'+' '+u'$\\mathrm{Ascension}$' + ' ' +u'$\\mathrm{(J2000)}$')
        ax1.axis_labels.set_ytext(u'$\\mathrm{Declination}$'+' '+u'$\\mathrm{(J2000)}$')
        ax1.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')
        ax1.set_tick_color('white')
        ax1.set_tick_labels_style('latex')
        ax1.set_labels_latex(True)
        ax1.show_rectangles(140.54131,51.002304,0.009504002,0.009504002, edgecolor='green', facecolor='none', lw=1)
        #ax1.show_beam(edgecolor='white', facecolor='none', lw=1)
        ax1.show_colorbar()
        ax1.colorbar.set_axis_label_text(u'$\\mathrm{Jy/beam}$')
        
        #Plot the spectrum
        ax2 = fig.add_axes([0.2,0.1,0.74,0.3])
        ax2.plot(velocity, spectrum)
        ax2.axvline(x=velocity[i],linewidth=2, color='green')
        ax2.axhline(y=0,color='grey')
        ax2.set_xlabel(u'$\\mathrm{Velocity (km/s)}$')
        ax2.set_ylabel(u'$\\mathrm{Flux}$'+u' '+u'$\\mathrm{Density (mJy/beam)}$')
        ax2.axis([max(velocity),min(velocity),-0.5,3.5])
        
        #Draw the figure and save as a PNG
	fig.canvas.draw()
	plt.savefig('pngs/'+str(i).zfill(3)+'.png', pad_inches=0)

#Create the animated GIF
os.system('convert -delay 10 -loop 100   pngs/*.png   animate.gif')
