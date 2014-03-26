#!python
"""
functions used for postprocessing my resonator characterisation data in the following shape:
columns of those text files:
up to case 205: f,|Z|,phase(I),hydro,mic,I,U
from case 206: f, U, I, phase(I), mic, phase(mic), hydro, phase(hydro), needle, phase(needle), |Z|, wavesamplemarker

Markus Stokmaier, KIT, IKET, March 2014
"""

from os import getcwd, listdir
from os.path import join, dirname, getmtime
from cPickle import Pickler
from pylab import sin, cos, exp, pi, sqrt, ceil, where, allclose, mean, fabs, argmax, argmin, absolute, cm
from pylab import array, asfarray, arange, ones, zeros, zeros_like, loadtxt, savetxt, size, shape, amin, amax
import numpy as np
from scipy import c_
import scipy.optimize as spo
import scipy.signal as sps
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ColorConverter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from time import time, localtime
from datetime import date as dtd


def give_datestring(tstamp=None):
    if tstamp is not None:
        return dtd.fromtimestamp(tstamp).strftime("%d. %B %Y")
    else:
        return dtd.fromtimestamp(time()).strftime("%d. %B %Y")

def parse_labview_time(timestamp):
    a,b=timestamp.split('.')
    h,m,s=a.split(':')
    h=int(float(h))
    m=int(float(m))
    s=int(float(s))
    frac=float(b)/10**len(b)
    return 3600*h + 60*m + s + frac

def rotate_phase(data,phi,bds=[-180,180]):
    """data is an array full of phase angles; phi will be added;
    the interval, given by bds, could also be [0,360] or [-pi,pi] or [-5*pi,-3*pi]"""
    data=asfarray(data)
    data+=phi
    data=cycle_into_domain(data,bds=bds)
    return data

def cycle_into_domain(data,bds=[-180,180]):
    lb,ub=asfarray(bds); assert lb<ub
    data=asfarray(data); w=ub-lb
    if any(data-lb<0):
        dist=ceil(where(data-lb<0,lb-data,0)/w)
        data+=dist*w
    if any(data-ub>0):
        dist=ceil(where(data-ub>0,data-ub,0)/w)
        data-=dist*w
    return data

def loglims(dat,lolim=None,hilim=None):
    """return something that makes sense for ax.set_ylim(list) in case of semilogy"""
    mn,mx=amin(dat),amax(dat)
    if mn<0:
        raise NotImplementedError()
    elif mn==0:
        if lolim is None:
            raise ValueError('please supply a wish for the lower limit, data goes down to zero')
        else:
            loval=lolim
    else:
        e1=int(float('{:e}'.format(mn)[-3:])); loval=10**e1
        e2=int(float('{:e}'.format(mx)[-3:]))+1; hival=10**e2
    return loval,hival

def deg2rad(d):
    return 2*pi*d/360.

def rad2deg(r):
    return 360.*r / (2*pi)

def expand_axis_for_pcolor(x,equidistant=False):
    """If you have two lists x and y of lengths n and m and an array dat of size (n,m) then pylab's pcolor
    will not display the last column and the last row. This little utility function helps you to create
    x2 and y2 so that they are one element longer and the pcolor patches will be centered."""
    if equidistant:
        n=len(x); d=(x[-1]-x[0])/float(n-1)
        return x[0]-0.5*d+arange(n+1,dtype=float)*d
    else:
        x2=[ x[0] - (x[1]-x[0])/2. ]
        for i in range(len(x)-1):
            x2.append((x[i]+x[i+1])/2.)
        x2.append( x[-1] + (x[-1]-x[-2])/2. )
        return array(x2)

def make_patch_spines_invisible(ax):
    """see http://matplotlib.sourceforge.net/examples/pylab_examples/multiple_yaxis_with_spines.html"""
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)




class DataContainer(object):

    def __init__(self,case,suffix='',lsmooth=True,datapath=None,switchiphase=False,phs_offset=None):
        self.case=case
        self.suffix=suffix
        loc=dirname(getcwd())
        nstr=str(case).zfill(2)
        self.hydroshow='Pascal'   # 'Pascal', 'bar', 'hyd'
        self.lsmooth=lsmooth
        self.switchiphase=switchiphase
        self.phs_offset=phs_offset
        prefix='case'
        if datapath is None:
            self.datapath=join(loc,'data')
        else:
            self.datapath=datapath
        self.data={} # for the data rows read from files
        self.edat={} # for electric data computed based upon U, I, and phase
        self.deciph={'voltage':'U', 'current':'I', 'phase':'I_phs', 'mic':'mic', 'hydro':'hyd', 'z':'Zamp'}
        self.read_data(prefix+nstr,suffix)
        self.process_basic()
        self.compute_gains()
        self.plotpath=join(loc,'plots')
        self.picklepath=join(loc,'pickles')
        self.textpath=join(loc,'analysis')
        self.outsuffix=''  # will be added to each plot and text file name (only output)
        self.axtitles=True
        self.dpi=120
        return

    def read_data(self,snippet,suffix):
        if suffix is None:
            fnp=join(self.datapath,snippet+'.txt')
        else:
            fnp=join(self.datapath,snippet+suffix+'.txt')
        rawdat=loadtxt(fnp)
        self.f=array(rawdat[:,0],copy=1)
        if 75 <= self.case < 206:
            dli=['Zamp','I_phs','hyd','mic','I','U']
        elif self.case >= 206:
            dli=['U','I','I_phs','mic','mic_phs','hyd','hyd_phs','ndl','ndl_phs','Zamp']
        else:
            raise ValueError('this subroutine should only be called for the modern files after case 74')
        for i,d in enumerate(dli):
            if self.switchiphase and d=='I_phs':
                self.data[d]=array(-rawdat[:,i+1],copy=1)
            else:
                self.data[d]=array(rawdat[:,i+1],copy=1)
            if (self.phs_offset is not None) and ('phs' in d):
                self.data[d]+=self.phs_offset
        self.ctime=getmtime(fnp)

    def data_cutting(self,finter,path,newsuffix):
        """
        save a subset defined by the frequency interval finter of the data
        """
        fini,fend=finter
        mark1=where(self.f>=fini,1,0)
        mark2=where(self.f<=fend,1,0)
        mark=mark1*mark2
        self.write_selected_data(mark,path,newsuffix)

    def write_selected_data(self,selected,writepath,newsuffix):
        """
        read in the raw data again and store it again in the writepath folder after having dumped all data
        points where the list "selected" contained a 0 or False.
        """
        if self.suffix is None:
            fname='case'+str(self.case).zfill(2)
        else:
            fname='case'+str(self.case).zfill(2)+self.suffix
        fnp=join(self.datapath,fname+'.txt')
        rawdat=loadtxt(fnp)
        assert shape(rawdat)[0]==len(selected)
        newdat=[]
        for i in xrange(len(selected)):
            if selected[i]:
                newdat.append(rawdat[i])
        savetxt(join(writepath,fname+newsuffix+'.txt'),array(newdat))

    def process_basic(self):
        if 'I_phs' in self.data: self.data['phsrad']=2*pi*self.data['I_phs']/360.                # phase in radian
        if 'U' in self.data: self.data['power']=self.data['U']*self.data['I']*cos(self.data['phsrad'])       # electric power input into piezo ring
        if 'hyd' in self.data: self.compute_bar_pressure()
        if ('U' in self.data) and ('I' in self.data) and ('I_phs' in self.data): self.compute_electric_data()

    def compute_gains(self):
        for sig in ['mic','hyd','ndl','pres','bpres']:
            if sig in self.data:
                self.data[sig[0]+'gain']=self.data[sig]/self.data['U']
                self.data[sig[0]+'gain']=where(self.data[sig[0]+'gain']!=np.inf,self.data[sig[0]+'gain'],-np.inf)

    def luxury_smooth(self):
        for sig in ['mic','hyd','ndl']:
            self.data[sig+'_sphs']=luxurious_smooth_finish(self.data[sig+'_phs'],nwin=32,offsteps=5,bds=[-180.,180.],order=4,win='hamming')

    def compute_bar_pressure(self):
        self.data['prespp']=self.data['hyd'] * 1000/10.06182 * 6894.75728 # 99.39 psi/V  --> 685239.6 Pa/V
        self.data['pres']=self.data['prespp']/2 # absolute peak (all values before that step were peak-to-peak)
        self.data['bpres']=self.data['pres']/1e5 # pressure in bar

    def tellme_soundpressure(self):
        i=argmax(self.data['hyd'])
        if self.data['bpres'][i] >= 0.1:
            r = r'max sound pressure is {:.3f} bar ({:.2f} psi), found at $\nu$ = {:.3f} kHz'
            r = r.format(self.data['bpres'][i], self.data['hyd'][i], self.f[i]/1e3)
        else:
            r = r'max sound pressure is {:.1f} Pa ({:.3f} psi), found at $\nu$ = {:.3f} kHz'
            r = r.format(self.data['pres'][i], self.data['hyd'][i], self.f[i]/1e3)
        return r

    def compute_electric_data(self):
        n=len(self.data['U'])
        u=zeros(n,dtype=complex)
        u.real[:]=self.data['U']
        i=cos(self.data['phsrad'])*self.data['I']+sin(self.data['phsrad'])*1j*self.data['I']
        z=u/i; y=1./z
        for var,lab in zip([u,i,z,y],['U','I','Z','Y']):
            self.edat[lab]=var

    def export_electric_data(self,path=None):
        data=zeros((len(self.f),9))
        data[:,0]=self.f
        for i,lab in enumerate(['U','I','Z','Y']):
            data[:,2*i+1]=self.edat[lab].real
            data[:,2*i+2]=self.edat[lab].imag
        hd='f Ure Uim Ire Iim R X G B'
        if path is None:
            fpn=join(self.picklepath,'edat_c'+str(self.case)+self.suffix+self.outsuffix+'.txt')
        else:
            fpn=join(path,'edat_c'+str(self.case)+self.suffix+self.outsuffix+'.txt')
        savetxt(fpn,data,header=hd,comments='')

    def standard_plot(self,ymax=None,finter=None):
        fig=plt.figure(figsize=[10,10])
        fig.subplots_adjust(right=0.82)
        ax1a=fig.add_subplot(311);             ax1b=ax1a.twinx(); ax1c=ax1a.twinx()
        if self.axtitles: plt.title('electric properties')
        ax2a=fig.add_subplot(312,sharex=ax1a); ax2b=ax2a.twinx(); ax2c=ax2a.twinx()
        if self.axtitles: plt.title('acoustic raw signals')
        ax3a=fig.add_subplot(313,sharex=ax1a); ax3b=ax3a.twinx(); ax3c=ax3a.twinx()
        if self.axtitles: plt.title('acoustic gain signals')
        make_patch_spines_invisible(ax1c)
        make_patch_spines_invisible(ax2c)
        make_patch_spines_invisible(ax3c)
        ax1c.spines["right"].set_position(("axes", 1.15))
        ax2c.spines["right"].set_position(("axes", 1.15))
        ax3c.spines["right"].set_position(("axes", 1.15))
        ax1c.spines["right"].set_visible(True)
        ax2c.spines["right"].set_visible(True)
        ax3c.spines["right"].set_visible(True)
        self.fill_standard_axes([ax1a, ax1b, ax1c, ax2a, ax2b, ax2c, ax3a, ax3b, ax3c], ymax=ymax, flim=finter)
        ax3a.set_xlabel(r'frequency $\nu$ in Hz')
        if self.axtitles:
            dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        if 'hyd' in self.data: fig.suptitle(self.tellme_soundpressure(),x=0.04,y=0.02,ha='left',va='bottom',fontsize=8)
        plt.savefig(join(self.plotpath,'standard_plot_c'+str(self.case).zfill(3)+self.suffix+'_D'+self.outsuffix+'.png'))
        plt.close()

    def fill_standard_axes(self,axlist,bright=1,flim=None,ymax=None):
        xlim_done=0; cols=[]; cli=['k','r','g','k','b','c','k','b','c']
        if self.hydroshow=='Pascal':
            dli=['U','I','I_phs','pres','mic','ndl','pgain','mgain','ngain']
            lli=['voltage $U$ in V','current $I$ in A','current phase $\phi(I)$ in degrees',
                 'hydrophone sound pressure in Pa','pill mic voltage','needle raw signal',
                 'hydrophone gain','pill mic gain','needle gain']
        elif self.hydroshow=='bar':
            dli=['U','I','I_phs','bpres','mic','ndl','bgain','mgain','ngain']
            lli=['voltage $U$ in V','current $I$ in A','current phase $\phi(I)$ in degrees',
                 'hydrophone sound pressure in barl','pill mic voltage','needle raw signal',
                 'hydrophone gain','pill mic gain','needle gain']
        elif self.hydroshow=='hyd':
            dli=['U','I','I_phs','hyd','mic','ndl','hgain','mgain','ngain']
            lli=['voltage $U$ in V','current $I$ in A','current phase $\phi(I)$ in degrees',
                 'hydrophone raw signal','pill mic voltage','needle raw signal',
                 'hydrophone gain','pill mic gain','needle gain']
        for i,d in enumerate(dli):
            if d in self.data:
                c1=brighten(ColConv.to_rgb(cli[i]),bright)
                axlist[i].plot(self.f, self.data[d], color=c1, ls='-', lw=0.7, alpha=0.5)
                cols.append(axlist[i].plot(self.f, self.data[d], mfc=c1, marker='o', ls='None', mec=c1, markersize=3, label=lli[i]))
                if 'phs' in d:
                    axlist[i].set_ylim(-180,180)
                    axlist[i].set_yticks([-180,-135,-90,-45,0,45,90,135,180])
                else:
                    axlist[i].set_ylim(bottom=0)
                axlist[i].set_ylabel(lli[i])
                axlist[i].yaxis.label.set_color(cli[i])
                if xlim_done==0:
                    if flim is not None:
                        axlist[i].set_xlim(flim)
                    else:
                        axlist[i].set_xlim([np.min(self.f),np.max(self.f)])
            else:
                cols.append(None)
                axlist[i].set_ylabel(lli[i])
                axlist[i].yaxis.label.set_color(cm.Greys(0.5))
        if ymax is not None:
            for ax,ym in zip(axlist,ymax):
                if ym is not None: ax.set_ylim(top=ym)

    def compare_with(self,othersdc,ymax=None,finter=None):
        fig=plt.figure(figsize=[10,10])
        fig.subplots_adjust(right=0.82)
        ax1a=fig.add_subplot(311);             ax1b=ax1a.twinx(); ax1c=ax1a.twinx()
        plt.title('electric properties')
        ax2a=fig.add_subplot(312,sharex=ax1a); ax2b=ax2a.twinx(); ax2c=ax2a.twinx()
        plt.title('acoustic raw signals')
        ax3a=fig.add_subplot(313,sharex=ax1a); ax3b=ax3a.twinx(); ax3c=ax3a.twinx()
        plt.title('acoustic gain signals')
        make_patch_spines_invisible(ax1c)
        make_patch_spines_invisible(ax2c)
        make_patch_spines_invisible(ax3c)
        ax1c.spines["right"].set_position(("axes", 1.15))
        ax2c.spines["right"].set_position(("axes", 1.15))
        ax3c.spines["right"].set_position(("axes", 1.15))
        ax1c.spines["right"].set_visible(True)
        ax2c.spines["right"].set_visible(True)
        ax3c.spines["right"].set_visible(True)
        if finter is None:
            fmin,fmax = min(self.f[0],othersdc.f[0]), max(self.f[-1],othersdc.f[-1])
        else:
            fmin,fmax=finter
        if self.hydroshow=='Pascal':
            dli=['U','I','I_phs','pres','mic','ndl','pgain','mgain','ngain']
        elif self.hydroshow=='bar':
            dli=['U','I','I_phs','bpres','mic','ndl','bgain','mgain','ngain']
        elif self.hydroshow=='hyd':
            dli=['U','I','I_phs','hyd','mic','ndl','hgain','mgain','ngain']
        if ymax is None:
            ymax=11*[0]
            for i,d in enumerate(dli):
                if 'phs' in d:
                    ymax[i]=None
                elif d not in self.data:
                    if d not in othersdc.data:
                        ymax[i]=None
                    else:
                        ymax[i]=np.max(othersdc.data[d])
                elif d not in othersdc.data:
                    ymax[i]=np.max(self.data[d])
                else:
                    ymax[i] = max(np.max(self.data[d]),np.max(othersdc.data[d]))
        #print 'ymax: ',ymax
        othersdc.fill_standard_axes([ax1a, ax1b, ax1c, ax2a, ax2b, ax2c, ax3a, ax3b, ax3c],2.5,flim=[fmin,fmax],ymax=ymax)
        self.fill_standard_axes([ax1a, ax1b, ax1c, ax2a, ax2b, ax2c, ax3a, ax3b, ax3c],flim=[fmin,fmax],ymax=ymax)
        ax3a.set_xlabel(r'frequency $\nu$ in Hz')
        odstr=dtd.fromtimestamp(othersdc.ctime).strftime("%d. %B %Y")
        dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
        txt='comparing case {} with {} (shaded), data from {} and {}'.format(self.case,othersdc.case,dstr,odstr)
        fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        txt = 'c '+str(self.case)+': '+self.tellme_soundpressure()+'\n' if 'hyd' in self.data else ''
        if 'hyd' in othersdc.data: txt+='c '+str(othersdc.case)+': '+othersdc.tellme_soundpressure()
        fig.suptitle(txt,x=0.04,y=0.02,ha='left',va='bottom',fontsize=8)
        plt.savefig(join(self.plotpath,'comp_c'+str(self.case).zfill(3)+self.suffix+'_with_c'+str(othersdc.case).zfill(3)+othersdc.suffix+'_D'+self.outsuffix+'.png'))
        plt.close()

    def electric_rawdat_plot(self,ymax=None,finter=None,logU=False,logI=False):
        fig=plt.figure(figsize=[8,4])
        ax1a=fig.add_subplot(111); ax1b=ax1a.twinx(); ax1c=ax1a.twinx()
        if self.axtitles: plt.title('raw data of electric measurements')
        make_patch_spines_invisible(ax1c)
        ax1c.spines["right"].set_position(("axes", 1.14))
        ax1c.spines["right"].set_visible(True)
        self.fill_electric_rawdat([ax1a, ax1b, ax1c], ymax=ymax, flim=finter)
        ax1a.set_xlabel(r'frequency $\nu$ in Hz')
        if logU:
            ax1a.semilogy(); ax1a.set_ylim(loglims(self.data['U']))
        if logI:
            ax1b.semilogy(); ax1b.set_ylim(loglims(self.data['I']))
        if self.axtitles:
            dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        plt.tight_layout()
        fig.subplots_adjust(right=0.75)
        plt.savefig(join(self.plotpath,'elrawdat_plot_c'+str(self.case).zfill(3)+self.suffix+self.outsuffix+'.png'),dpi=self.dpi)
        plt.close()

    def fill_electric_rawdat(self,axlist,bright=1,flim=None,ymax=None):
        xlim_done=0; cols=[]; cli=['k','r','g']
        dli=['U','I','I_phs']
        lli=['voltage $U$ in V','current $I$ in A','current phase $\phi(I)$ in degrees']
        for i,d in enumerate(dli):
            if d in self.data:
                c1=brighten(ColConv.to_rgb(cli[i]),bright)
                axlist[i].plot(self.f, self.data[d], color=c1, ls='-', lw=0.7, alpha=0.5)
                cols.append(axlist[i].plot(self.f, self.data[d], mfc=c1, marker='o', ls='None', mec=c1, markersize=3, label=lli[i]))
                if 'phs' in d:
                    axlist[i].set_ylim(-180,180)
                    axlist[i].set_yticks([-180,-135,-90,-45,0,45,90,135,180])
                else:
                    axlist[i].set_ylim(bottom=0)
                axlist[i].set_ylabel(lli[i])
                axlist[i].yaxis.label.set_color(cli[i])
                if xlim_done==0:
                    if flim is not None:
                        axlist[i].set_xlim(flim); #print 'setting xlim: ',flim
                    else:
                        axlist[i].set_xlim([np.min(self.f),np.max(self.f)])
            else:
                cols.append(None)
                axlist[i].set_ylabel(lli[i])
                axlist[i].yaxis.label.set_color(cm.Greys(0.5))
        if ymax is not None:
            for ax,ym in zip(axlist,ymax):
                if ym is not None: ax.set_ylim(top=ym)

    def compare_electric_rawdat_plot(self,othersdc,ymax=None,finter=None,logU=False,logI=False):
        fig=plt.figure(figsize=[8,4])
        fig.subplots_adjust(right=0.75)
        ax1a=fig.add_subplot(111); ax1b=ax1a.twinx(); ax1c=ax1a.twinx()
        #if self.axtitles: plt.title('raw data of electric measurements')
        make_patch_spines_invisible(ax1c)
        ax1c.spines["right"].set_position(("axes", 1.14))
        ax1c.spines["right"].set_visible(True)
        othersdc.fill_electric_rawdat([ax1a, ax1b, ax1c],2.5, ymax=ymax, flim=finter)
        self.fill_electric_rawdat([ax1a, ax1b, ax1c], ymax=ymax, flim=finter)
        ax1a.set_xlabel(r'frequency $\nu$ in Hz')
        if logU:
            ax1a.semilogy(); ax1a.set_ylim(loglims(self.data['U']))
        if logI:
            ax1b.semilogy(); ax1b.set_ylim(loglims(self.data['I']))
        if self.axtitles:
            #dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            #odstr=dtd.fromtimestamp(othersdc.ctime).strftime("%d. %B %Y")
            #txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            txt='comparing case {} with {} (shaded)'.format(self.case,othersdc.case)
            #txt='comparing case {} with {} (shaded), data from {} and {}'.format(self.case,othersdc.case,dstr,odstr)
            fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        #plt.tight_layout()
        plt.savefig(join(self.plotpath,'elrawdat_comp_c'+str(self.case).zfill(3)+self.suffix+'_with_c'+str(othersdc.case).zfill(3)+othersdc.suffix+self.outsuffix+'.png'))
        plt.close()

    def electric_plot(self,ymin=None,ymax=None,equalize=False,flipup=False):
        fig=plt.figure(figsize=[10,10])
        ax1a=fig.add_subplot(311);             ax1b=ax1a.twinx()
        if self.axtitles: plt.title(r'magnitudes of impedance $|Z|=|R+iX|$ and admittance $|Y|=|G+iB|$')
        ax2a=fig.add_subplot(312,sharex=ax1a); ax2b=ax2a.twinx()
        if self.axtitles: plt.title(r'resistance $R$ and reactance $X$')
        ax3a=fig.add_subplot(313,sharex=ax1a); ax3b=ax3a.twinx()
        if self.axtitles: plt.title('conductance $G$ and susceptance $B$')
        self.fill_electric_axes([ax1a, ax1b, ax2a, ax2b, ax3a, ax3b], ymin=ymin, ymax=ymax, equalize=equalize, flipup=flipup)
        ax3a.set_xlabel(r'frequency $\nu$ in Hz')
        if self.axtitles:
            dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        plt.savefig(join(self.plotpath,'electric_plot_c'+str(self.case).zfill(3)+self.suffix+self.outsuffix+'.png'),dpi=self.dpi)
        plt.close()

    def compare_power(self,othersdc,ymin=None,ymax=None,logy=False):
        fig=plt.figure(figsize=[8,4])
        ax=fig.add_subplot(111)
        grey=brighten(ColConv.to_rgb('k'),2.5)
        #ax.plot(othersdc.f,othersdc.data['power'],fc=grey,ec=grey,marker='o')
        ax.plot(othersdc.f, othersdc.data['power'], mfc=grey, marker='o', ls='None', mec=grey, markersize=3)
        ax.plot(self.f,self.data['power'],'ko')
        #if self.axtitles: plt.title('dissipated power')
        ax.set_xlabel(r'frequency $\nu$ in Hz')
        ax.set_ylabel(r'dissipated power $P$ in Watt')
        if logy:
            ax.semilogy()
            if (ymin is None) and (ymax is None):
                ax.set_ylim(loglims([min(amin(self.data['power']),amin(othersdc.data['power'])),max(amax(self.data['power']),amax(othersdc.data['power']))]))
        if ymin is not None: ax.set_ylim(ymin=ymin)
        if ymax is not None: ax.set_ylim(ymax=ymax)
        if self.axtitles:
            dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        plt.tight_layout()
        plt.savefig(join(self.plotpath,'comp_power_c'+str(self.case).zfill(3)+self.suffix+'_with_c'+str(othersdc.case).zfill(3)+othersdc.suffix+self.outsuffix+'.png'))
        plt.close()

    def fill_electric_axes(self,axlist,bright=1,flim=None,ymin=None,ymax=None,equalize=False,flipup=False):
        cols=[]; cli=['c','r','b','purple','orange','crimson']
        if flipup:
            dli=[absolute(self.edat['Z']), absolute(self.edat['Y']),
                 fabs(self.edat['Z'].real), fabs(self.edat['Z'].imag),
                 fabs(self.edat['Y'].real), fabs(self.edat['Y'].imag)]
        else:
            dli=[absolute(self.edat['Z']), absolute(self.edat['Y']),
                 self.edat['Z'].real, self.edat['Z'].imag,
                 self.edat['Y'].real, self.edat['Y'].imag]
        lli=['$|Z|$ in $\Omega$','$|Y|$ in $1/\Omega$','$R$ in $\Omega$','$X$ in $\Omega$','$G$ in $1/\Omega$','$B$ in $1/\Omega$']
        #print 'ymin: ',ymin
        #print 'ymax: ',ymax
        for i,d in enumerate(dli):
            c1=brighten(ColConv.to_rgb(cli[i]),bright)
            axlist[i].plot(self.f, d, color=c1, ls='-', lw=0.7, alpha=0.5)
            cols.append(axlist[i].plot(self.f, d, mfc=c1, marker='o', ls='None', mec=c1, markersize=3, label=lli[i]))
            if flipup: axlist[i].set_ylim(bottom=0)
            axlist[i].set_ylabel(lli[i])
            axlist[i].yaxis.label.set_color(cli[i])
            if flim is not None:
                axlist[i].set_xlim(flim)
            else:
                axlist[i].set_xlim([np.min(self.f),np.max(self.f)])
            if ymax is not None:
                if ymax[i] is not None:
                    axlist[i].set_ylim(top=ymax[i])
                    #print 'pushing ymax={} to {}'.format(ymax[i],axlist[i])
            if ymin is not None:
                if ymin[i] is not None:
                    axlist[i].set_ylim(bottom=ymin[i])
                    #print 'pushing ymin={} to {}'.format(ymin[i],axlist[i])
        if equalize:
            equalize_ymax([axlist[0],axlist[2],axlist[3]])
            equalize_ymax([axlist[1],axlist[4],axlist[5]])
            if not flipup:
                equalize_ymin([axlist[0],axlist[2],axlist[3]])
                equalize_ymin([axlist[1],axlist[4],axlist[5]])

    def compare_electric(self,othersdc,finter=None,ymin=None,ymax=None,equalize=False,flipup=False):
        fig=plt.figure(figsize=[10,10])
        fig.subplots_adjust(right=0.75)
        ax1a=fig.add_subplot(311);             ax1b=ax1a.twinx()
        plt.title('electric properties')
        ax2a=fig.add_subplot(312,sharex=ax1a); ax2b=ax2a.twinx()
        plt.title('acoustic raw signals')
        ax3a=fig.add_subplot(313,sharex=ax1a); ax3b=ax3a.twinx()
        plt.title('acoustic gain signals')
        #if finter is None:
        #    fmin,fmax = min(self.f[0],othersdc.f[0]), max(self.f[-1],othersdc.f[-1])
        #    fini=0; fend=-1; ofini=0; ofend=-1
        #else:
        #    fmin,fmax=finter
        #    fini=int(argmin(fabs(self.f-finter[0]))); fend=int(argmin(fabs(self.f-finter[1])))
        #    ofini=int(argmin(fabs(othersdc.f-finter[0]))); ofend=int(argmin(fabs(othersdc.f-finter[1])))
        #fini,fend=finifend(self.f,finter)
        #ofini,ofend=finifend(othersdc.f,finter)
        #print 'fini,fend, ofini,ofend: ',fini,fend,ofini,ofend
        if ymax is None:
            ymax=6*[0]
            ymax[0] = relevant_max_duo(self.f,absolute(self.edat['Z']),othersdc.f,absolute(othersdc.edat['Z']),finter)
            ymax[1] = relevant_max_duo(self.f,absolute(self.edat['Y']),othersdc.f,absolute(othersdc.edat['Y']),finter)
            ymax[2] = relevant_max_duo(self.f,self.edat['Z'].real,othersdc.f,othersdc.edat['Z'].real,finter)
            ymax[3] = relevant_max_duo(self.f,self.edat['Z'].imag,othersdc.f,othersdc.edat['Z'].imag,finter)
            ymax[4] = relevant_max_duo(self.f,self.edat['Y'].real,othersdc.f,othersdc.edat['Y'].real,finter)
            ymax[5] = relevant_max_duo(self.f,self.edat['Y'].imag,othersdc.f,othersdc.edat['Y'].imag,finter)
            #ymax[0] = max(np.nanmax(absolute(self.edat['Z'][fini:fend])),np.nanmax(absolute(othersdc.edat['Z'][ofini:ofend])))
            #ymax[1] = max(np.nanmax(absolute(self.edat['Y'][fini:fend])),np.nanmax(absolute(othersdc.edat['Y'][ofini:ofend])))
            #ymax[2] = max(np.nanmax(self.edat['Z'][fini:fend].real),np.nanmax(othersdc.edat['Z'][ofini:ofend].real))
            #ymax[3] = max(np.nanmax(self.edat['Z'][fini:fend].imag),np.nanmax(othersdc.edat['Z'][ofini:ofend].imag))
            #ymax[4] = max(np.nanmax(self.edat['Y'][fini:fend].real),np.nanmax(othersdc.edat['Y'][ofini:ofend].real))
            #ymax[5] = max(np.nanmax(self.edat['Y'][fini:fend].imag),np.nanmax(othersdc.edat['Y'][ofini:ofend].imag))
        if ymin is None:
            ymin=6*[0]
            if not flipup:
                ymin[0] = relevant_min_duo(self.f,absolute(self.edat['Z']),othersdc.f,absolute(othersdc.edat['Z']),finter)
                ymin[1] = relevant_min_duo(self.f,absolute(self.edat['Y']),othersdc.f,absolute(othersdc.edat['Y']),finter)
                ymin[2] = relevant_min_duo(self.f,self.edat['Z'].real,othersdc.f,othersdc.edat['Z'].real,finter)
                ymin[3] = relevant_min_duo(self.f,self.edat['Z'].imag,othersdc.f,othersdc.edat['Z'].imag,finter)
                ymin[4] = relevant_min_duo(self.f,self.edat['Y'].real,othersdc.f,othersdc.edat['Y'].real,finter)
                ymin[5] = relevant_min_duo(self.f,self.edat['Y'].imag,othersdc.f,othersdc.edat['Y'].imag,finter)
                #ymin[0] = min(np.nanmin(absolute(self.edat['Z'][fini:fend])),np.nanmin(absolute(othersdc.edat['Z'][ofini:ofend])))
                #ymin[1] = min(np.nanmin(absolute(self.edat['Y'][fini:fend])),np.nanmin(absolute(othersdc.edat['Y'][ofini:ofend])))
                #ymin[2] = min(np.nanmin(self.edat['Z'][fini:fend].real),np.nanmin(othersdc.edat['Z'][ofini:ofend].real))
                #ymin[3] = min(np.nanmin(self.edat['Z'][fini:fend].imag),np.nanmin(othersdc.edat['Z'][ofini:ofend].imag))
                #ymin[4] = min(np.nanmin(self.edat['Y'][fini:fend].real),np.nanmin(othersdc.edat['Y'][ofini:ofend].real))
                #ymin[5] = min(np.nanmin(self.edat['Y'][fini:fend].imag),np.nanmin(othersdc.edat['Y'][ofini:ofend].imag))
        if finter is None:
            fmin,fmax = min(self.f[0],othersdc.f[0]), max(self.f[-1],othersdc.f[-1])
        else:
            fmin,fmax=finter
        othersdc.fill_electric_axes([ax1a, ax1b, ax2a, ax2b, ax3a, ax3b],2.5, flim=[fmin,fmax] ,ymax=ymax, ymin=ymin)
        self.fill_electric_axes([ax1a, ax1b, ax2a, ax2b, ax3a, ax3b], flim=[fmin,fmax] ,ymax=ymax, ymin=ymin)
        if equalize:
            equalize_ymax([ax1a,ax2a,ax2b])
            equalize_ymax([ax1b,ax3a,ax3b])
            if not flipup:
                equalize_ymin([ax1a,ax2a,ax2b])
                equalize_ymin([ax1b,ax3a,ax3b])
        ax3a.set_xlabel(r'frequency $\nu$ in Hz')
        odstr=dtd.fromtimestamp(othersdc.ctime).strftime("%d. %B %Y")
        dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
        txt='comparing case {} with {} (shaded), data from {} and {}'.format(self.case,othersdc.case,dstr,odstr)
        fig.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        txt = 'c '+str(self.case)+': '+self.tellme_soundpressure()+'\n' if 'hyd' in self.data else ''
        if 'hyd' in othersdc.data: txt+='c '+str(othersdc.case)+': '+othersdc.tellme_soundpressure()
        fig.suptitle(txt,x=0.04,y=0.02,ha='left',va='bottom',fontsize=8)
        #plt.savefig(join(self.plotpath,'ref_c'+str(othersdc.case).zfill(3)+'_comp_c'+str(self.case).zfill(3)+'_B.png'),dpi=self.dpi)
        plt.savefig(join(self.plotpath,'electric_comp_c'+str(self.case).zfill(3)+self.suffix+'_with_c'+str(othersdc.case).zfill(3)+othersdc.suffix+self.outsuffix+'.png'),dpi=self.dpi)
        plt.close()

    def real_i_per_u_plot(self,ymin=None,ymax=None,logy=True):
        plt.plot(self.f,self.data['I']/self.data['U'],'kx',label='raw data I/U')
        plt.plot(self.f,1./self.data['Zamp'],'y-',label='raw data 1/|Z|')
        plt.plot(self.f,self.edat['I'].real/self.edat['U'].real,'cd',label='el. data Re(I)/Re(U)')
        if logy: plt.semilogy()
        if ymin is not None: plt.ylim(ymin=ymin)
        if ymax is not None: plt.ylim(ymax=ymax)
        plt.xlabel('frequency in Hz')
        plt.ylabel('Re(I)/Re(U) in Siemens')
        plt.legend()
        if self.axtitles:
            dstr=dtd.fromtimestamp(self.ctime).strftime("%d. %B %Y")
            txt='case {}, chamber characterisation data from {}'.format(self.case,dstr)
            plt.suptitle(txt,x=0.5,y=0.98,ha='center',va='top',fontsize=12)
        plt.savefig(join(self.plotpath,'real_i_per_u_plot_c'+str(self.case).zfill(3)+self.suffix+self.outsuffix+'.png'),dpi=self.dpi)
        plt.close()

    def small_printout(self):
        ofl=open(join(self.textpath,'printout_c'+str(self.case).zfill(2)+self.suffix+self.outsuffix+'.txt'),'w')
        for key,dat in sorted(self.data.iteritems()):
            idx1=np.argmin(dat); idx2=np.argmax(dat)
            txt='{:12s} ranges from {:>30} to {:>30}'.format(key,dat[idx1],dat[idx2])
            txt+=' with min at f= {:>16} and max at f= {:>16}\n'.format(self.f[idx1],self.f[idx2])
            ofl.write(txt)
        ofl.close()

    def determine_Q(self,lab,finter,nwin_pd,nwin_s=8,smoothing='standard',qmode='sqrt'):
        """
        arguments: lab -> which data label, e.g. 'mic'; finter -> frequency interval [a,b];
        nwin -> length of window for smoothing; nwin_pd -> window for peak detection.
        Attention, function peak_q_factors() from below is used here which doesn't expect
        any noise. This means for you that you should make sure by looking at the checkplot, that the output makes sense.
        qmode -> mode of Q-factor calculation, if 'half' then FWHM, if 'sqrt' then where decayed to 1/sqrt(2) of max value
        (qmode='sqrt' should be used for all amplitudes which give something proportional to power when squared, 'half' should
        be used if the signal itself is proportional to power)
        """
        i1=argmin(fabs(self.f-finter[0])); i2=argmin(fabs(self.f-finter[1]))
        x=self.f[i1:i2]; y=self.data[lab][i1:i2]
        if smoothing=='standard':
            b, a = sps.butter(4, 1/(len(y)/20.))
            sy = sps.filtfilt(b, a, y)
        elif smoothing=='luxury':
            sy=luxurious_smooth_finish(y,nwin=nwin_s)
        elif smoothing=='none':
            sy=y
        else:
            raise ValueError("argument smoothing must be in ['standard', 'luxury', 'none']")
        qidx,q=peak_q_factors(x, sy, nwin_pd=nwin_pd, mode=qmode); #print qidx; print q
        ofl=open(join(self.textpath,'Q_factors_c'+str(self.case).zfill(2)+self.suffix+'_'+lab+'.txt'),'w')
        ofl.write('scanned frequency interval: {} to {} Hz\n\n'.format(x[0],x[-1]))
        for ii,qq in zip(qidx,q):
            if type(qq)==str:
                ofl.write('peak at {} Hz yielded error: "{}" (max value = {})\n'.format(x[ii],qq,sy[ii]))
            else:
                ofl.write('peak at {} Hz has Q = {} (max value = {})\n'.format(x[ii],qq,sy[ii]))
        ofl.close()
        """ Now let's still make the plot. """
        plt.plot(x,y,'yo')
        plt.plot(x,sy,'b-',lw=2)
        plt.title('case {}{}: Q-factors of peaks in {} data'.format(self.case,self.suffix,lab))
        txt=''
        for ii,qq in zip(qidx,q):
            if type(qq)==str:
                plt.axhline(y=sy[ii],color='r')
                plt.axvline(x=x[ii],color='r')
                plt.annotate(qq, xy=(x[ii],sy[ii]),  xycoords='data', xytext=(60, -40), textcoords='offset points', size=14, #bbox=dict(boxstyle="round", fc="0.8"),
                             arrowprops=dict(arrowstyle="simple", fc='b', ec="none", connectionstyle="arc3,rad=0.3", shrinkB=15, alpha=0.4))
                txt+='peak at f = {:.1f} Hz yielded error: "{}" (max value = {:.3e})\n'.format(x[ii],qq,sy[ii])
            else:
                plt.axhline(y=sy[ii],color='g')
                plt.axvline(x=x[ii],color='g')
                atxt=r'$Q={:.3f}$'.format(qq)
                plt.annotate(atxt, xy=(x[ii],sy[ii]),  xycoords='data', xytext=(60, -40), textcoords='offset points', size=14, #bbox=dict(boxstyle="round", fc="0.8"),
                             arrowprops=dict(arrowstyle="simple", fc='b', ec="none", connectionstyle="arc3,rad=0.3", shrinkB=15, alpha=0.4))
                txt+='peak at f = {:.1f} has a Q-factor of {:.3f} (max value = {:.3e})\n'.format(x[ii],qq,sy[ii])
        plt.suptitle(txt[:-2],x=0.03,y=0.01,ha='left',va='bottom',fontsize=8)
        plt.savefig(join(self.plotpath,'Q_factors_c'+str(self.case).zfill(2)+self.suffix+'_'+lab+self.outsuffix+'.png'),dpi=self.dpi)
        plt.close()

    def freq_of_max(self,varname,finter=None):
        if varname not in self.data:
            raise ValueError("choose a variable in SDC_with_dict.data")
        if finter is None:
            fidx=argmax(self.data[varname])
        else:
            i1,i2=[argmin(fabs(self.f-f)) for f in finter]
            fidx=argmax(self.data[varname][i1:i2])+i1
        return fidx,self.f[fidx]
        

    def pickle_self(self):
        ofile=open(join(self.picklepath,'sdc_pickle_c'+str(self.case).zfill(2)+self.suffix+self.outsuffix+'.txt'), 'w')
        einmachglas=Pickler(ofile)
        einmachglas.dump(self)
        ofile.close()



def argrelmax1d(data, nwin=1):
    """for 1D-sequences"""
    N=len(data)
    results = np.ones(N, dtype=bool)
    for shift in xrange(1, nwin + 1):
        plus = np.roll(data,shift)
        minus = np.roll(data,-shift)
        results &= np.greater(data, plus)
        results &= np.greater(data, minus)
    idx=[]
    for i,val in enumerate(results[nwin:-nwin]):
        if val: idx.append(i+nwin)
    return idx

def peak_q_factors(f,a,nwin_pd=25,mode='sqrt'):
    """
    finds peaks and determines Q-factors on smooth curves with no noise
    keyword argument nwin_pd: window size for relative maximum determination
    returns: list of peak indices, list of Q-factors
    returns error string inside the Q-factor list if
    a) peak has a side peak with a little minimum in between
    b) peak's asymmetry exceeds some threshold
    possible improvement to implement:
    shoulder as third sort of error to be detected
    """
    f=np.asfarray(f); a=np.asfarray(a); N=len(f); assert len(a)==N
    rm=argrelmax1d(a,nwin=nwin_pd)
    qidx=[]; q=[]
    for idx in rm:
        maxval=a[idx]; val=maxval; ls,rs=0,0
        if mode == 'half':
            hval=0.5*maxval
        elif mode == 'sqrt':
            hval=maxval/sqrt(2.)
        upwardsbend=False
        sidecrash=False
        while val > hval:
            ls+=1
            if idx-ls < 0:
                sidecrash=True
                break
            val=a[idx-ls]
            if a[idx-ls] >= a[idx-ls+1]:
                upwardsbend=True
                break # if it goes upwards again
        val=maxval
        while val > hval:
            rs+=1
            if idx+rs > N-1:
                sidecrash=True
                break
            val=a[idx+rs]
            if a[idx+rs] >= a[idx+rs-1]:
                upwardsbend=True
                break # if it goes upwards again
        if sidecrash:
            qidx.append(idx); q.append('sidecrash')
        elif upwardsbend:
            qidx.append(idx); q.append('upbend')
        else:
            fres=f[idx]
            fleft=f[idx-ls] + ((hval-a[idx-ls])/(a[idx-ls+1]-a[idx-ls])) * (f[idx-ls+1]-f[idx-ls])
            fright=f[idx+rs] + ((hval-a[idx+rs])/(a[idx+rs-1]-a[idx+rs])) * (f[idx+rs-1]-f[idx+rs])
            assert fleft<fres
            assert fright>fres
            if not 1./2.5 < (fres-fleft)/(fright-fres) < 2.5:
                qidx.append(idx); q.append('asym')
            else:
                qidx.append(idx); q.append(fres/(fright-fleft))
    return qidx,q


def equalize_ymax(axlist):
    ymax=[]
    for ax in axlist:
        ymax.append(ax.get_ylim()[1])
    themax=amax(ymax); #print 'maxima for ymax equaliser: ',ymax
    for ax in axlist:
        ax.set_ylim(top=themax)

def equalize_ymin(axlist):
    ymin=[]
    for ax in axlist:
        ymin.append(ax.get_ylim()[0])
    themin=amin(ymin); #print 'minima for ymin equaliser: ',ymin
    for ax in axlist:
        ax.set_ylim(bottom=themin)

def interval_touched(x,i,min_overlap=2):
    """check whether (some part of) the interval i lies within the range covered by x;
    x is a monotonically ascending sequence of values; i=[a,b] with a<=b is the interval"""
    mo=min_overlap
    if mo < len(x): return False
    elif (x[-mo] < i[0]) or (x[mo] > i[-1]): return False
    else: return True

def separate_intervals(x,i):
    if (x[-1]<i[0]) or (x[0]>i[-1]): return True
    else: return False

def dat_iv_overlap(data,interval):
    assert len(interval)==2
    idx1=int(argmin(fabs(data-interval[0])))
    idx2=int(argmin(fabs(data-interval[1])))
    if idx1==idx2:
        if (data[idx1]<interval[0]) or (data[idx1]>interval[1]):
            return 0,[]
        else:
            return 1,[idx1]
    else:
        #return (idx2-idx1) * (int(idx2>idx1)*2-1) # |idx2-idx1|
        assert idx2>idx1
        return idx2-idx1+1, [idx1,idx2]


def relevant_min_uno(x,y,xinter):
    if xinter is None:
        return np.nanmin(y)
    else:
        ovl,lims=dat_iv_overlap(x,xinter)
        if ovl==0:
            return None
        elif ovl==1:
            return y[lims[0]]
        else:
            return np.nanmin(y[lims[0]:lims[1]])

def relevant_max_uno(x,y,xinter):
    if xinter is None:
        return np.nanmax(y)
    else:
        ovl,lims=dat_iv_overlap(x,xinter)
        if ovl==0:
            return None
        elif ovl==1:
            return y[lims[0]]
        else:
            return np.nanmax(y[lims[0]:lims[1]])

def relevant_min_duo(x1,y1,x2,y2,xinter):
    mn1=relevant_min_uno(x1,y1,xinter)
    mn2=relevant_min_uno(x2,y2,xinter)
    if (mn1 is None) and (mn2 is None):
        return None
    elif (mn1 is None):
        return mn2
    elif (mn2 is None):
        return mn1
    else:
        return min(mn1,mn2)

def relevant_max_duo(x1,y1,x2,y2,xinter):
    mx1=relevant_max_uno(x1,y1,xinter)
    mx2=relevant_max_uno(x2,y2,xinter)
    if (mx1 is None) and (mx2 is None):
        return None
    elif (mx1 is None):
        return mx2
    elif (mx2 is None):
        return mx1
    else:
        return max(mx1,mx2)



#--- polynomial fitting utilities

class ParabolicFitTools:
    def __init__(self,xdata,ydata):
        self.xdata=xdata
        self.ydata=ydata
        self.n=len(xdata)
    def model(self,params):
        a,b,c=params
        return a + b*self.xdata + c*self.xdata**2
    def rescalc(self,params):
        return self.ydata-self.model(params)
    def jacobian(self,params):
        return c_[-ones(self.n),-self.xdata,-self.xdata**2]

class PolyFitTools:
    def __init__(self,xdata,ydata,order=2):
        self.xdata=xdata
        self.ydata=ydata
        self.n=len(xdata)
        self.o=order
    def model(self,params):
        if len(params)!=self.o+1:
            msg='this PolyFitTools instance of fitting order {0} does not suit the call with a length {1} vector'.format(self.o,len(params))
            raise TypeError(msg)
        result=zeros(self.n)
        for i,par in enumerate(params):
            result+=par*self.xdata**i
        return result
    def rescalc(self,params):
        return self.ydata-self.model(params)
    def jacobian(self,params):
        Df=-ones((self.n,self.o+1))
        for i in range(1,self.o+1):
            Df[:,i]*=self.xdata**i
        return Df

def fit_data(xdata,ydata,order=2):
    pft=PolyFitTools(xdata,ydata,order=order)
    xstart=array([mean(ydata)]+order*[0])
    bestpars,info=spo.leastsq(pft.rescalc,xstart,Dfun=pft.jacobian)
    result=zeros_like(ydata)
    for i,par in enumerate(bestpars):
        result+=par*xdata**i
    return result

def offset_fit(xdata,ydata,offset=0.,bds=[-180.,180.],order=2):
    y=ydata+offset
    y=cycle_into_domain(ydata+offset,bds=bds)
    fit=fit_data(xdata,y,order=order)
    return cycle_into_domain(fit-offset,bds=bds)

#def closeness_in_cyclic_domain(x,y,bds=[-180.,180.]):
#    bw=max(bds)-min(bds); dist=zeros((3,len(y)))
#    for i,shift in enumerate([-bw,0.,bw]):
#        dist[i,:]=(y+shift-x)**2
#    mindist=np.min(dist,axis=0)
#    return np.sum(mindist)

def closeness_in_cyclic_domain(x,y,bds=[-180.,180.]):
    bw=max(bds)-min(bds)
    if len(shape(x))==1:
        dist=zeros((3,len(y)))
        for i,shift in enumerate([-bw,0.,bw]):
            dist[i,:]=(y+shift-x)**2
        mindist=np.min(dist,axis=0)
        return np.sum(mindist)
    elif len(shape(x))==2:
        n1,n2=shape(x); dist=zeros((3,n1,n2))
        for i,shift in enumerate([-bw,0.,bw]):
            dist[i,:,:]=(y+shift-x)**2
        mindist=np.min(dist,axis=0)
        return np.sum(mindist)
    else:
        raise NotImplementedError('this function can only deal with 1D or 2D arrays')

def RMS_distance_in_cyclic_domain(x,y,bds=[-180.,180.]):
    bw=max(bds)-min(bds)
    if len(shape(x))==1:
        sqdist=zeros((3,len(y)))
        for i,shift in enumerate([-bw,0.,bw]):
            sqdist[i,:]=(y+shift-x)**2
        minsqdist=np.min(sqdist,axis=0)
        return sqrt(np.sum(minsqdist)/size(x))
    elif len(shape(x))==2:
        n1,n2=shape(x); sqdist=zeros((3,n1,n2))
        for i,shift in enumerate([-bw,0.,bw]):
            sqdist[i,:,:]=(y+shift-x)**2
        minsqdist=np.min(sqdist,axis=0)
        return sqrt(np.sum(minsqdist)/size(x))
    else:
        raise NotImplementedError('this function can only deal with 1D or 2D arrays')

def multi_offset_fit(xdata,ydata,steps=3,bds=[-180.,180.],addoffset=0.,order=2):
    offset=(max(bds)-min(bds))/float(steps)
    fitlist=[]; qualities=[]
    for i in range(steps):
        fitlist.append(offset_fit(xdata,ydata,offset=i*offset+addoffset,order=order))
        qualities.append(closeness_in_cyclic_domain(fitlist[-1],ydata,bds=bds))
    best=array(qualities).argmin()
    return fitlist[best]

def multi_offset_mean(data,steps=3,bds=[-180.,180.],addoffset=0.):
    offset=(max(bds)-min(bds))/float(steps); n=len(data)
    meanlist=[]; qualities=[]
    for i in range(steps):
        shifted_data=cycle_into_domain(data+i*offset+addoffset)
        meanlist.append(np.mean(shifted_data)-i*offset+addoffset)
        qualities.append(closeness_in_cyclic_domain(meanlist[-1]*ones(n),data,bds=bds))
    meanlist=cycle_into_domain(array(meanlist),bds=bds)
    best=array(qualities).argmin()
    return meanlist[best]

def luxurious_smooth_finish(data,nwin=8,offsteps=3,bds=[-180.,180.],addoffset=0.,order=2,win='flat'):  #,sidestep=1):
    """
    A smoothing filter for a 1D data set based on qubic fits of pieces. To be used when the signal
    values cover a cyclic domain as happening with phase data.
    The filter is quite costly, but it can handle one particular nastiness often happening
    when treating the phase data of an oscillating signal: when the noisy phase signal leaves a bounded
    domain (e.g. from 0 to 2*pi or from -180 to 180 degrees) and comes in again from the other side it
    sometimes hesitantly jumps back and forth a couple times. Where other filters will smear the signal
    so it hovers a lot around the center of the domain, this filter will conserve/produce a sharp jump
    to the other side because it is built on the assumption of a cyclic domain of values. Make sure the
    data window size (argument nwin) is small enough so the data that can fall within it never traverses
    the domain multiple times.
    The keyword argument offsteps determines the number of offset steps each signal piece is tried with,
    thus calculation cost is proportional to it. The offset is used to cycle data through the domain, so
    doing this a couple more times in smaller steps increases the probability that a data segment is
    away from the domain boundaries at least once. Hence, 3, 4, and 5 are probably the most reasonable
    settings, although for signals known to stay away from part of the domain or to move slowly with low
    noise using less offset steps will completely suffice (in combination with having chosen an
    appropriate additional constant offset in the case of offsteps=1).
    """
    N=len(data); n=nwin; fitlib=zeros((n,N)); fit=zeros(N); x=arange(n)
    for i in range(N-n+1):
        y=data[i:i+n]
        fitpiece=multi_offset_fit(x,y,steps=offsteps,bds=bds,addoffset=addoffset,order=order)
        for j in range(n):
            fitlib[n-j-1,i+j]=fitpiece[j]
    if win in [None,'None','flat']:
        for i in range(N):
            iidx=max(n-i-1,0); fidx=min(N-i,n)
            fit[i]=multi_offset_mean(fitlib[iidx:fidx,i],steps=offsteps,bds=bds,addoffset=addoffset)
    else:
        for i in range(n-1):
            # first we must do something against the triangles full of zeros at each end of fitlib
            fitlib[:n-1-i,i]=multi_offset_mean(fitlib[n-1-i:,i],steps=offsteps,bds=bds,addoffset=addoffset)
            fitlib[i+1:,-i-1]=multi_offset_mean(fitlib[:i+1,-i-1],steps=offsteps,bds=bds,addoffset=addoffset)
        bw=max(bds)-min(bds); shifts=[-bw,0.,bw]; dist=zeros((n,3)); choices=zeros(n,dtype=int)
        if win=='bartlett': w=np.bartlett(nwin)
        elif win=='blackman': w=np.blackman(nwin)
        elif win=='hamming': w=np.hamming(nwin)
        elif win=='hanning': w=np.hanning(nwin)
        elif win=='flat': w=ones(nwin)
        else: raise TypeError("the keyword argument 'win' allows only those entries: None, 'None', 'flat', 'bartlett', 'blackman', 'hamming', 'hanning'")
        w/=mean(w)
        for i in range(N):
            fit[i]=multi_offset_mean(fitlib[:,i],steps=offsteps,bds=bds,addoffset=addoffset)
            dist[:,0]=fabs(fitlib[:,i]+shifts[0]-fit[i]); dist[:,1]=fabs(fitlib[:,i]+shifts[1]-fit[i]); dist[:,2]=fabs(fitlib[:,i]+shifts[2]-fit[i])
            choices[:]=[dist[j,:].argmin() for j in range(n)]
            fitlib[:,i]+=[shifts[choices[j]] for j in range(n)]
            fitlib[:,i]*=w
        fit[:]=cycle_into_domain(mean(fitlib,axis=0),bds=bds)
    return fit

def roll_towards_continuity(data,bds=[-180.,180.],largest_delta_fraction=0.2):
    """
    Imagine phase data jumping from -180 to 180 degrees, i.e. leaving the domain boudaries through the bottom and reentering
    the region from above. This function will roll down the data so the signal does not care any more about the boundaries
    and just continuously continues its way down. The signal should be rather smooth and noise low. The argument bds is a list
    containing lower and upper boundary; largest_delta_fraction is the largest expected jump from one data point to the next
    divided by the width of the whole domain.
    """
    N=len(data); ldf=largest_delta_fraction; bw=max(bds)-min(bds)
    diff=data[1:]-data[:-1]
    downcrossings=where(diff>(1-ldf)*bw,1,0)
    upcrossings=where(diff<(ldf-1)*bw,1,0)
    crossings=upcrossings-downcrossings
    shift=bw*asfarray([0]+[np.sum(crossings[:i+1]) for i in range(N-1)])
    return data+shift

#--- custom colormaps
#--- a) blue_red
# a color map from blue over almost neutrally yellowish white to red
cdict4 = {'red':  ((0.0, 0.0, 0.0),(0.25,0.0, 0.0),(0.5, 1.0, 1.0),(0.75,1.0, 1.0),(1.0, 0.3, 1.0)),
         'green': ((0.0, 0.0, 0.0),(0.25,0.0, 0.0),(0.5, 1.0, 1.0),(0.75,0.0, 0.0),(1.0, 0.0, 0.0)),
         'blue':  ((0.0, 0.0, 0.3),(0.25,1.0, 1.0),(0.5, 0.5, 0.5),(0.75,0.0, 0.0),(1.0, 0.0, 0.0))}
blue_red4 = LinearSegmentedColormap('BlueRed4', cdict4)

cdict5a = {'red':  ((0.0, 0.0, 0.0),(0.25,0.0, 0.0),(0.5, 1.0, 1.0),(0.75,1.0, 1.0),(1.0, 0.3, 1.0)),
          'green': ((0.0, 0.0, 0.2),(0.25,0.4, 0.0),(0.5, 1.0, 1.0),(0.75,0.0, 0.3),(1.0, 0.1, 0.0)),
          'blue':  ((0.0, 0.0, 0.3),(0.25,1.0, 1.0),(0.5, 0.5, 0.5),(0.75,0.0, 0.0),(1.0, 0.0, 0.0))}
blue_red5a = LinearSegmentedColormap('BlueRed5a', cdict5a)

cdict5b = {'red':  ((0.0, 0.0, 0.0),(0.25,0.0, 0.0),(0.5, 1.0, 1.0),(0.75,1.0, 1.0),(1.0, 0.3, 1.0)),
          'green': ((0.0, 0.0, 0.0),(0.25,0.0, 0.4),(0.5, 1.0, 1.0),(0.75,0.3, 0.0),(1.0, 0.0, 0.0)),
          'blue':  ((0.0, 0.0, 0.3),(0.25,1.0, 1.0),(0.5, 0.5, 0.5),(0.75,0.0, 0.0),(1.0, 0.0, 0.0))}
blue_red5b = LinearSegmentedColormap('BlueRed5b', cdict5b)

cdict5c = {'red':  ((0.0, 0.0, 0.0),(0.25,0.0, 0.0),(0.5, 1.0, 1.0),(0.75,1.0, 1.0),(1.0, 0.3, 1.0)),
          'green': ((0.0, 0.0, 0.0),(0.25,0.0, 0.2),(0.5, 1.0, 1.0),(0.75,0.2, 0.0),(1.0, 0.0, 0.0)),
          'blue':  ((0.0, 0.0, 0.3),(0.25,1.0, 1.0),(0.5, 0.5, 0.5),(0.75,0.0, 0.0),(1.0, 0.0, 0.0))}
blue_red5c = LinearSegmentedColormap('BlueRed5c', cdict5c)

cdict6b = {'red':  ((0.0, 0.0, 0.0),(0.125,0.00, 0.00),(0.25,0.0, 0.1),(0.375,0.00, 0.00),(0.5, 1.0, 1.0),(0.625,1.00, 1.00),(0.75,0.8, 1.0),(0.875,0.65, 0.75),(1.0, 0.3, 0.0)),
          'green': ((0.0, 0.0, 0.0),(0.125,0.20, 0.00),(0.25,0.0, 0.0),(0.375,0.75, 0.40),(0.5, 1.0, 1.0),(0.625,0.40, 0.75),(0.75,0.0, 0.0),(0.875,0.00, 0.20),(1.0, 0.0, 0.0)),
          'blue':  ((0.0, 0.0, 0.3),(0.125,0.75, 0.65),(0.25,1.0, 0.8),(0.375,0.75, 0.80),(0.5, 1.0, 1.0),(0.625,0.20, 0.25),(0.75,0.1, 0.0),(0.875,0.00, 0.00),(1.0, 0.0, 0.0))}
blue_red6b = LinearSegmentedColormap('BlueRed6b', cdict6b)

cdict6d = {'red':  ((0.0, 0.0, 0.0),(0.125,0.00, 0.00),(0.25,0.3, 0.0),(0.375,0.00, 0.00),(0.5, 1.0, 1.0),(0.625,1.00, 1.00),(0.75,0.8, 1.0),(0.875,0.35, 0.75),(1.0, 0.2, 0.0)),
          'green': ((0.0, 0.0, 0.0),(0.125,0.30, 0.00),(0.25,0.0, 0.4),(0.375,0.75, 0.40),(0.5, 1.0, 1.0),(0.625,0.40, 0.75),(0.75,0.2, 0.0),(0.875,0.00, 0.30),(1.0, 0.0, 0.0)),
          'blue':  ((0.0, 0.0, 0.2),(0.125,0.75, 0.35),(0.25,1.0, 0.4),(0.375,0.75, 0.80),(0.5, 1.0, 1.0),(0.625,0.20, 0.25),(0.75,0.0, 0.1),(0.875,0.00, 0.00),(1.0, 0.0, 0.0))}
blue_red6d = LinearSegmentedColormap('BlueRed6d', cdict6d)

cdict7 = {'red':  ((0.0,0.3, 0.3),(0.5,0.2, 0.2),(1.0, 0.8, 0.0)),
         'green': ((0.0,0.1, 0.1),(0.5,0.2, 0.2),(1.0, 1.0, 0.0)),
         'blue':  ((0.0,0.0, 0.0),(0.5,1.0, 1.0),(1.0, 0.9, 0.0))}
blue_red7 = LinearSegmentedColormap('BlueRed7', cdict7)

#--- b) sidekick
cdict_sidekick = {'red':  ((0.00, 1.0, 1.0),(0.12, 0.9, 0.9),(0.25, 0.5, 0.5),(0.33, 0.1, 0.1),(0.40, 0.2, 0.2),(0.45, 0.4, 0.4),(0.50, 1.0, 1.0),(0.55, 0.0, 0.0),(0.60, 0.1, 0.1),(0.67, 0.1, 0.1),(0.75, 0.0, 0.0),(0.86, 0.0, 0.0),(1.00, 0.9, 0.9)),
                  'green':((0.00, 1.0, 1.0),(0.12, 0.4, 0.4),(0.25, 0.1, 0.1),(0.33, 0.0, 0.0),(0.40, 0.0, 0.0),(0.45, 0.0, 0.0),(0.50, 1.0, 1.0),(0.55, 0.4, 0.4),(0.60, 0.2, 0.2),(0.67, 0.1, 0.1),(0.75, 0.4, 0.4),(0.86, 0.7, 0.7),(1.00, 1.0, 1.0)),
                  'blue': ((0.00, 0.7, 0.7),(0.12, 0.0, 0.0),(0.25, 0.0, 0.0),(0.33, 0.0, 0.0),(0.40, 0.3, 0.3),(0.45, 0.0, 0.0),(0.50, 1.0, 1.0),(0.55, 0.4, 0.4),(0.60, 0.2, 0.2),(0.67, 0.5, 0.5),(0.75, 0.2, 0.2),(0.86, 0.0, 0.0),(1.00, 0.2, 0.2))}
sidekick = LinearSegmentedColormap('sidekick', cdict_sidekick)

#--- custom color stuff

ColConv=ColorConverter()

def brighten(rgb,factor):
    rgb=asfarray(rgb)
    if factor > 1:
        wdist=1-rgb    # how far is each channel from being full
        wdist/=factor  # if factor is 3, then bring the color 3 times closer to white
        return tuple(1-wdist)
    elif 0 < factor < 1:
        rgb*=factor  # if factor is 0.3, then reduce each channel to 30%
        return tuple(rgb)
    elif factor==1:
        return tuple(rgb)
    else:
        raise ValueError("You cant't go into the negative range with the rgb channels.")