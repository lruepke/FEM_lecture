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
import meshio
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LightSource
from nice import niceAxis,text3d

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
        ax.text(x[el].mean(),-0.15,'Element %d'%(i),ha='center',va='top')
        # node
        ax.scatter(x[el],x[el]*0-0.1,marker='o',fc=colors_node,ec='w',s=100,clip_on=False,zorder=10)
        ax.text(x[el][0],-0.1, '%d'%(el[0]),ha='center',va='center',zorder=11,fontsize=9,color='w')
        if(i==(len(el2nod)-1)): 
            ax.text(x[el][1],-0.1, '%d'%(el[1]),ha='center',va='center',zorder=11,fontsize=9,color='w')
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

def loadMesh2D(gmshfile):
    mesh=meshio.read(gmshfile)
    GCOORD=mesh.points[:,0:2]
    el2node=[]
    for cells in mesh.cells:
        if(cells.type=='quad'):
            el2node=cells.data
    return {'GCOORD':GCOORD, 'EL2NOD':el2node}
def plotMesh_2D(mesh):
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    el2nod=mesh['EL2NOD']
    len_x,len_y=x.max()-x.min(),y.max()-y.min()
    figwidth=12
    figheight=len_y/len_x*figwidth
    fig=plt.figure(figsize=(figwidth,figheight))
    ax=plt.gca()
    # plot element
    for element in el2nod:
        ind=np.append(element,element[0])
        ax.fill(x[ind],y[ind])
    # plot node
    for nod in np.unique(el2nod.reshape(-1)):
        ax.plot(x[nod],y[nod],marker='o',ms=20,mfc='w',mec='k')
        ax.text(x[nod],y[nod],'%d'%(nod+1),ha='center',va='center')
    plt.show()
# shape function
def shapes(s1, s2):
    N1 = 0.25*(1-s1)*(1-s2)
    N2 = 0.25*(1+s1)*(1-s2)
    N3 = 0.25*(1+s1)*(1+s2)
    N4 = 0.25*(1-s1)*(1+s2)
    return N1, N2, N3, N4
def write2VTU(vtufile, xx,yy,zz):
    ncols=xx.shape[0]
    nrows=xx.shape[1]
    x=xx.reshape(-1)
    y=yy.reshape(-1)
    z=zz.reshape(-1)
    npoints=len(x)
    nCells=(ncols-1)*(nrows-1)
    VTK_CELLTYPE=9 #四边形
    np_per_cell=4
    fpout=open(vtufile,'w')
    fpout.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
    fpout.write('  <UnstructuredGrid>\n')
    fpout.write('    <Piece NumberOfPoints="%.0f" NumberOfCells="%.0f">\n'%(npoints,nCells))
    fpout.write('      <PointData Scalars="ShapeFunction">\n')
    fpout.write('        <DataArray type="Float64" Name="ShapeFunction" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing ShapeFunction field ...')
    for i in range(0,len(z)):
        fpout.write('%f '%(z[i]))
    fpout.write('\n        </DataArray>\n')
    fpout.write('      </PointData>\n')
    fpout.write('      <CellData>\n')
    fpout.write('      </CellData>\n')
    fpout.write('      <Points>\n')
    fpout.write('        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">\n')
    # Info('Writing xyz data ...')
    for i in range(0,len(x)):
        fpout.write('          %f %f %f\n'% (x[i],y[i],z[i]))
    fpout.write('        </DataArray>\n')
    fpout.write('      </Points>\n')
    fpout.write('      <Cells>\n')
    fpout.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
    # Info('Writing connectivity ...')
    for nrow in range(0,nrows-1):
        for ncol in range(0,ncols-1):
            LL=ncol + nrow*ncols
            fpout.write('          %.0f %.0f %.0f %.0f\n'%(LL, LL+1, LL+1+ncols, LL+ncols))
    fpout.write('        </DataArray>\n')
    fpout.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing offsets ...')
    for i in range(0,nCells):
        fpout.write('%.0f '%(i*np_per_cell))
    fpout.write('        </DataArray>\n')
    fpout.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing cell type ...')
    for i in range(0,nCells):
        fpout.write('%.0f '%(VTK_CELLTYPE))
    fpout.write('        </DataArray>\n')
    fpout.write('      </Cells>\n')
    fpout.write('    </Piece>\n')
    fpout.write('  </UnstructuredGrid>\n')
    fpout.write('</VTKFile>\n')
    # Info('xyz to vtu Done')
    fpout.close()
    # Info('Converting ASCII to binary')
    os.system('meshio-binary '+vtufile)
def Linear2D_Quad(xmin=0,dx=1,ymin=0,dy=1):
    nip     = 4
    gauss   = np.array([[-1, 1, 1, -1], [-1, -1, 1, 1]]) * np.sqrt(1.0/3.0)
    shapes(gauss[0,1],gauss[1,1])
    x=np.linspace(-1,1,50)
    y=np.linspace(-1,1,50)
    xx,yy=np.meshgrid(x,y)
    N1, N2, N3, N4=shapes(xx,yy)
    # fig,axes=plt.subplots(2,2)
    # for ax,N in zip([axes[0][0],axes[0][1],axes[1][0],axes[1][1]], [N1, N2, N3, N4]):
    #     CS=ax.contourf(xx,yy,N, cmap='jet')
    #     plt.colorbar(CS)
    # plt.show()

    # write2VTU('N1.vtu',xx,yy,N1)
    # write2VTU('N2.vtu',xx,yy,N2)
    # write2VTU('N3.vtu',xx,yy,N3)
    # write2VTU('N4.vtu',xx,yy,N4)
    sign_shapefunc=[['-','-'],['+','-'],['+','+'],['-','+']]
    formula_shapefunc = lambda signs : '$\\frac{1}{4}(1%s\\xi)(1%s\\eta)$'%(signs[0],signs[1])
    fig = plt.figure(figsize=(16,8))
    for i, N in enumerate([N1, N2, N3, N4]):
        ax = fig.add_subplot(1,4,i+1, projection='3d',facecolor='None')
        ls = LightSource(270, 45)
        rgb = ls.shade(N, cmap=cm.plasma, vert_exag=0.1, blend_mode='soft')
        CS=ax.plot_surface(xx, yy, N, rstride=1, cstride=1, facecolors=rgb,
        linewidth=0, antialiased=True, shade=False)
        ax.set_xlabel('$\\xi$',labelpad=0)
        ax.set_ylabel('$\\eta$',labelpad=0)
        ax.set_zlabel('Shape function',labelpad=-3)
        ax.zaxis.set_minor_locator(MultipleLocator(0.1))
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(0,1)
        # ax.xaxis.set_ticks([])
        # ax.yaxis.set_ticks([])
        ax.zaxis.set_ticks([0,1])
        ax.set_title('N$_{\mathregular{%d}}=$%s'%(i,formula_shapefunc(sign_shapefunc[i])))
        # 重新自定义坐标轴属性
        niceAxis(ax,fill_pane=False,label3D=True,fs_label=0.1, scaled=False)
    plt.subplots_adjust(wspace=0)
    for fmt in fmt_figs:
        figname=str('%s/shapeFunction_2D_Q1.%s'%(figpath,fmt))
        plt.savefig(figname, bbox_inches='tight')
def main(argv):
    # Linear1D()
    Linear2D_Quad()
if __name__ == '__main__':
    sys.exit(main(sys.argv))