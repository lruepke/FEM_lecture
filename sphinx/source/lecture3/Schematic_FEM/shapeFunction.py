#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Plot schematic diagram of shape functions of FEM
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/05/25, GEOMAR
# ===============================================================

import sys
import argparse
import sys
import os
from colored import fg, bg, attr
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')
#===============================================================
import linecache
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
figpath='.'
fmt_figs=['pdf','svg']

def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Plot schematic diagram of shape functions of FEM'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/05/25, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))
def createMesh_1D(xmin=0,xmax=10,dx=2):
    x=np.arange(xmin,xmax,dx)
    el2nod=[]
    for i in range(0,len(x)-1):
        el2nod.append([i, i+1])
    return {'GCOORD':{'x':x}, 'EL2NOD':el2nod}

def Linear1D():
    colors=mpl.cm.get_cmap('tab10').colors
    dx=2
    mesh=createMesh_1D(dx=2)
    x=mesh['GCOORD']['x']
    el2nod=mesh['EL2NOD']
    # figure
    fig=plt.figure(figsize=(8,2.5))
    ax=plt.gca()
    # shape function definition
    Ni_x  = lambda xi, dx, x : 1 - (x-xi)/dx
    Ni1_x = lambda xi, dx, x : (x-xi)/dx
    # plot mesh and shape function
    for i, el in enumerate(el2nod):
        colors_node=[colors[el[0]%len(colors)], colors[el[1]%len(colors)]]
        # element
        ax.plot(x[el],x[el]*0-0.1,lw=4,clip_on=False,color='k')
        ax.text(x[el].mean(),-0.15,'Element %d'%(i+1),ha='center',va='top')
        # node
        ax.scatter(x[el],x[el]*0-0.1,marker='o',fc=colors_node,ec='w',s=100,clip_on=False,zorder=10)
        ax.text(x[el][0],-0.1, '%d'%(el[0]+1),ha='center',va='center',zorder=11,fontsize=9,color='w')
        if(i==(len(el2nod)-1)): 
            ax.text(x[el][1],-0.1, '%d'%(el[1]+1),ha='center',va='center',zorder=11,fontsize=9,color='w')
        # calculate shape function value of each node
        tmp_x=np.linspace(x[el][0],x[el][1],2)
        Ni = Ni_x(x[el][0], dx, tmp_x)
        Ni1 = Ni1_x(x[el][0], dx, tmp_x)
        ax.plot(tmp_x, Ni,   color=colors_node[0])
        ax.plot(tmp_x, Ni1,  color=colors_node[1])
    # axis settings
    ax.set_ylim(0,1)
    ax.set_xlim(x.min(),x.max())
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.xaxis.set_ticks([])
    ax.spines['left'].set_position(("axes", -0.02))
    ax.set_ylabel('Shape function value')
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    for fmt in fmt_figs:
        figname=str('%s/shapeFunction_Linear_1D.%s'%(figpath,fmt))
        plt.savefig(figname, bbox_inches='tight')

def main(argv):
    Linear1D()

if __name__ == '__main__':
    sys.exit(main(sys.argv))