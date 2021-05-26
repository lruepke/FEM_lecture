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
def createMesh_2D(xmin=0,dx=1,nx=5,ymin=0,dy=1,ny=4):
    x=np.arange(xmin,xmin+(nx)*dx,dx)
    y=np.arange(ymin,ymin+(ny)*dy,dy)
    xx,yy=np.meshgrid(x,y)
    xx,yy=xx.reshape(-1,1),yy.reshape(-1,1)
    GCOORD=np.column_stack((xx,yy))
    # el2nod
    EL2NOD=[]
    for j in range(0,ny-1):
        for i in range(0,nx-1):
            LL=i+(nx)*j
            EL2NOD.append(np.array([LL, LL+1, LL+nx+1, LL+nx],dtype=int))
    EL2NOD=np.array(EL2NOD,dtype=int)
    # print(GCOORD.shape,EL2NOD)
    return {'GCOORD':GCOORD, 'EL2NOD':EL2NOD}
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
        figurename=str('%s/shapeFunction_Linear_1D.%s'%(figpath,fmt))
        plt.savefig(figurename, bbox_inches='tight')

def loadMesh2D(gmshfile):
    mesh=meshio.read(gmshfile)
    GCOORD=mesh.points[:,0:2]
    EL2NOD=[]
    for cells in mesh.cells:
        if((cells.type=='quad') | (cells.type=='triangle')):
            EL2NOD=cells.data
    return {'GCOORD':GCOORD, 'EL2NOD':EL2NOD}
def plotMesh_2D(mesh,colorsName='tab20',figname='mesh2D_structured',extend_x=0):
    colors=mpl.cm.get_cmap(colorsName).colors
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    el2nod=mesh['EL2NOD']
    len_x,len_y=x.max()-x.min(),y.max()-y.min()
    figwidth=12
    figheight=len_y/len_x*figwidth
    fig=plt.figure(figsize=(figwidth,figheight))
    ax=plt.gca()
    # plot element
    for i,element in enumerate(el2nod):
        ind=np.append(element,element[0])
        ax.fill(x[ind],y[ind],color=colors[i%len(colors)])
        ax.text(x[element].mean(),y[element].mean(),'Element %d'%(i),ha='center',va='center',bbox={'fc':'lightgray','ec':'None','boxstyle':'round'})
    # plot node
    nodes=[]
    for el in el2nod:
        for node in el:
            nodes.append(node)
    nodes=np.array(nodes,dtype=int)
    for nod in np.unique(nodes):
        ax.plot(x[nod],y[nod],marker='o',ms=15,mfc='w',mec='k',clip_on=False)
        ax.text(x[nod],y[nod],'%d'%(nod),ha='center',va='center')
    ax.axis('off')
    ax.axis('scaled')
    ax.set_xlim(x.min()-len_x*extend_x,x.max()+len_x*extend_x)
    ax.set_ylim(y.min(),y.max())
    for fmt in fmt_figs:
        figurename=str('%s/%s.%s'%(figpath,figname,fmt))
        plt.savefig(figurename, bbox_inches='tight')
def plotConnectivityMatrix_2D(mesh,colorsName='tab20',figname='Matrix2D_structured'):
    colors=mpl.cm.get_cmap(colorsName).colors
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    el2nod=mesh['EL2NOD']
    len_x,len_y=x.max()-x.min(),y.max()-y.min()
    nnel,nnod=len(el2nod),len(x)
    nnel1=int(np.sqrt(nnel))
    nnel2=int(nnel/nnel1)
    rows=np.min([nnel1,nnel2])
    cols=int(nnel/rows)
    width_ratios=[1]*(cols*2+1)
    width_ratios[cols]=0.1
    figheight=12
    figwidth=figheight*2.1
    fig,axes=plt.subplots(rows,cols*2+1, sharex=True, sharey=True,figsize=(figwidth,figheight),
    gridspec_kw={"hspace":0.15,'wspace':0.15,"width_ratios":width_ratios})
    # plot matrix of each element
    gs=axes[0][0].get_gridspec()
    ax_global = fig.add_subplot(gs[:, cols+1:])
    ax_equal = fig.add_subplot(gs[:,cols])
    for ax in axes[:,cols:]:
        for a in ax:
            a.remove()
    for i in range(0,rows):
        for j in range(0,cols):
            ind_el = j+i*cols
            element=el2nod[ind_el]
            ax_local=axes[i][j]
            for ax in [ax_local,ax_global]:
                ax.axis('scaled')
                ax.set_xlim(-0.5,nnod-0.5)
                ax.set_ylim(-0.5,nnod-0.5)
                ax.xaxis.set_ticks(np.arange(0,nnod+1,1)-0.5)
                ax.yaxis.set_ticks(np.arange(0,nnod+1,1)-0.5)
                ax.grid(axis='both',which='major',clip_on=False)
                ax.tick_params(axis='both',which='both',color='None')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                for spine in ax.spines:
                    ax.spines[spine].set_visible(False)
                ax.invert_yaxis()
                ax_local.set_title('Element %d'%(ind_el))
                for nod_col in element:
                    y=[nod_col-0.5, nod_col-0.5, nod_col+0.5, nod_col+0.5]
                    for nod_row in element:
                        x=[nod_row-0.5, nod_row+0.5, nod_row+0.5, nod_row-0.5]
                        ax.fill(x,y,color=colors[ind_el%len(colors)],alpha=0.8)
                if(j<(cols-1)):
                    ax_local.plot(1.06,0.5,'P',color='k',transform=ax_local.transAxes,clip_on=False)
                if((i>0) & (j==0)):
                    ax_local.plot(-0.06,0.5,'P',color='k',transform=ax_local.transAxes,clip_on=False)
    ax_equal.plot([0,1],[0.505,0.505],lw=3,color='k',transform=ax_equal.transAxes,clip_on=False)
    ax_equal.plot([0,1],[0.495,0.495],lw=3,color='k',transform=ax_equal.transAxes,clip_on=False)
    ax_equal.axis('off')
    # for i,element in enumerate(el2nod):
    #     ind=np.append(element,element[0])
    #     axes[0][0].fill(x[ind],y[ind],color=colors[i%len(colors)])
        # ax.text(x[element].mean(),y[element].mean(),'Element %d'%(i),ha='center',va='center',bbox={'fc':'lightgray','ec':'None','boxstyle':'round'})
    # axes[0].remove()
    for fmt in fmt_figs:
        figurename=str('%s/%s.%s'%(figpath,figname,fmt))
        plt.savefig(figurename, bbox_inches='tight')
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
        figurename=str('%s/shapeFunction_2D_Q1.%s'%(figpath,fmt))
        plt.savefig(figurename, bbox_inches='tight')
def connectivity():
    # 1. structured mesh with structured indexing
    mesh=createMesh_2D()
    plotMesh_2D(mesh,extend_x=1/3.0)
    plotConnectivityMatrix_2D(mesh)
    # structured mesh with unstructured indexing
    mesh=loadMesh2D('mesh/structure.msh')
    plotMesh_2D(mesh,figname='mesh2D_structured_usi',extend_x=1/3.0)
    plotConnectivityMatrix_2D(mesh,figname='Matrix2D_structured_usi')
    # unstructured mesh
    mesh=loadMesh2D('mesh/unstructure.msh')
    plotMesh_2D(mesh,figname='mesh2D_unstructured')
    plotConnectivityMatrix_2D(mesh,figname='Matrix2D_unstructured')
    
def main(argv):
    # Linear1D()
    Linear2D_Quad()
    # connectivity()
if __name__ == '__main__':
    sys.exit(main(sys.argv))