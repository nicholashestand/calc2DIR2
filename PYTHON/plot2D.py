#!/usr/bin/env python
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib

def getUserArgs():
    """Read command line arguments and return as an argparse type""" 
    parser = argparse.ArgumentParser(description='Plot 2D contour plot from file specified by user.')
    parser.add_argument('file', metavar='file', nargs='+', help='the name of the file containing the data')
    parser.add_argument('--columns', '-cols', metavar='x y z', type=int, nargs=3, help='which columns in file are x y and z data')
    parser.add_argument('--xbounds', '-xb', metavar='x0 x1', type=float, nargs=2, help='lower and upper x range')
    parser.add_argument('--ybounds', '-yb', metavar='y0 y1', type=float, nargs=2, help='lower and upper y range')
    parser.add_argument('--xlabel', '-xl', metavar='x label', type=str, nargs=1, help='x label')
    parser.add_argument('--ylabel', '-yl', metavar='y label', type=str, nargs=1, help='y label')
    parser.add_argument('--label', '-l', metavar='label', type=str, nargs='+', help='label' )
    parser.add_argument('--output', '-o', metavar='output file', type=str, nargs=1, help='output file name')
    return parser
 
args = getUserArgs()
fns_in = args.parse_args().file
cols = args.parse_args().columns if args.parse_args().columns else [0,1,2]

nplt=len(fns_in)
fig, ax = plt.subplots(nplt,figsize=(3.3,nplt*2.6),sharex=True,sharey=True,squeeze=False)

for i in range(nplt):
    fn_in = fns_in[i]
    print("Open file {}.".format( fn_in ))
    data=np.loadtxt(fn_in)
    print("Done reading file.")

    for col in cols:
        if col > len(data[0,:]):
            raise ValueError('Column {} is larger than the number of columns in {}.'.format(col,fn_in))

    print ("Plotting {}.".format(fn_in))
    # normalize data to max abs z value
    maxv=data[:,cols[2]].max()
    minv=data[:,cols[2]].min()
    maxa=max( abs(maxv), abs(minv) )
    data[:,cols[2]] /= maxa
    nlevels=18
    fontsize=8
    lw=0.5
    levels=np.linspace(1./nlevels,1,num=nlevels, endpoint=True )
    levels=np.concatenate((-1.*np.flipud(levels),levels))

    dimx=np.unique(data[:,cols[0]]).size
    dimy=np.unique(data[:,cols[1]]).size
    Xdata=data[:,cols[0]].reshape((dimx,dimy))
    Ydata=data[:,cols[1]].reshape((dimx,dimy))
    Zdata=data[:,cols[2]].reshape((dimx,dimy))

    ax[i,0].contourf( Xdata, Ydata, Zdata, levels, cmap=plt.get_cmap("bwr") )
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    ax[i,0].contour( Xdata, Ydata, Zdata, levels, colors='k', linewidths=lw )
    if ( args.parse_args().xbounds ):
        xlim = args.parse_args().xbounds
        ax[i,0].set_xlim(xlim)
    if ( args.parse_args().ybounds ):
        ylim = args.parse_args().ybounds
        ax[i,0].set_ylim(ylim)
    if ( args.parse_args().xlabel ):
        ax[i,0].set_xlabel(args.parse_args().xlabel[0], fontsize=fontsize)
    if ( args.parse_args().ylabel ):
        ax[i,0].set_ylabel(args.parse_args().ylabel[0], fontsize=fontsize)
    if ( args.parse_args().label ):
        a=(args.parse_args().label[i],)
        ax[i,0].legend(a,fontsize=fontsize, frameon=False)

    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='both', labelsize=fontsize, direction='in', top=True,bottom=True, left=True, right=True)

plt.tight_layout(pad=0.1)
if ( args.parse_args().output ):
    plt.savefig(args.parse_args().output[0])
else:
    plt.show()
