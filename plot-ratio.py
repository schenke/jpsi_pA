from math import *
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np

font = {'color':  'black',
    'weight': 'normal',
    'size': 18
    }
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Palatino']
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelpad'] = 15

with PdfPages('plots-jpsi.pdf') as pdf:

    fig = plt.figure(figsize=(5,3.8))
    ax = plt.gca()

    plt.title(r'Forward: p-going',weight='bold',fontsize=16)

    # generic setup
    ax.grid(True, lw=0.5, ls=':', zorder=0)
    plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )

    plt.xlabel(r'$Q_{s,p}$[GeV]',fontsize=16)

    plt.ylabel(r'$J/\Psi$ yield',fontsize=16)

    
    av_g, avg_error, avj_fwd, avj_fwd_error, avj_bwd, avj_bwd_error, av_qsp,av_qst,av_fwd_qsp, av_fwd_qst, av_bwd_qsp, av_bwd_qst  = np.loadtxt('averages_mv.dat', delimiter =' ', unpack = True)
    glue, glue_error, jpsi_fwd, jpsierror_fwd, jpsi_bwd, jpsierror_bwd, qsp, qst, qsp_fwd, qst_fwd,qsp_bwd, qst_bwd  = np.loadtxt('qsdep-mv.dat', delimiter =' ', unpack = True)
    
    glue_yield = glue/av_g
    jpsi_yield_fwd = jpsi_fwd/avj_fwd
    jpsi_yield_bwd = jpsi_bwd/avj_bwd
    
   # lims = np.linspace(0,15)
   # plt.plot(lims,lims,color='grey',ls='--',lw=1.5)
    plt.scatter(qsp_fwd,jpsi_yield_fwd,c=qst_fwd,s=40,marker='o', cmap='viridis')
    cbar = plt.colorbar()
    cbar.set_label('$Q_{s,A}$[GeV]', rotation=270,fontsize=16) 
    pdf.savefig(bbox_inches='tight') 
    plt.close()

    fig = plt.figure(figsize=(5,3.8))
    ax = plt.gca()

    plt.title(r'Backward: Pb-going',weight='bold',fontsize=16)

    # generic setup
    ax.grid(True, lw=0.5, ls=':', zorder=0)
    plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )

    plt.xlabel(r'$Q_{s,p}$[GeV]',fontsize=16)

    plt.ylabel(r'$J/\Psi$ yield',fontsize=16)
    
    jpsi_yield_bwd = jpsi_bwd/avj_bwd
    
   # lims = np.linspace(0,15)
   # plt.plot(lims,lims,color='grey',ls='--',lw=1.5)
    plt.scatter(qsp_bwd,jpsi_yield_bwd,c=qst_bwd,s=40,marker='o', cmap='viridis')
    cbar = plt.colorbar()
    cbar.set_label('$Q_{s,A}$[GeV]', rotation=270,fontsize=16) 
    pdf.savefig(bbox_inches='tight') 
    plt.close()

    fig = plt.figure(figsize=(5,3.8))
    ax = plt.gca()

    plt.title(r'Mid-rapidity',weight='bold',fontsize=16)

    # generic setup
    ax.grid(True, lw=0.5, ls=':', zorder=0)
    plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )

    plt.xlabel(r'$Q_{s,p}$[GeV]',fontsize=16)

    plt.ylabel(r'Gluon yield',fontsize=16)
    
    jpsi_yield_bwd = jpsi_bwd/avj_bwd
    
   # lims = np.linspace(0,15)
   # plt.plot(lims,lims,color='grey',ls='--',lw=1.5)
    plt.scatter(qsp,jpsi_yield_bwd,c=qst,s=40,marker='o', cmap='viridis')
    cbar = plt.colorbar()
    cbar.set_label('$Q_{s,A}$[GeV]', rotation=270,fontsize=16) 
    pdf.savefig(bbox_inches='tight') 
    plt.close()