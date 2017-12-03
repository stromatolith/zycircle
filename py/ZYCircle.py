#!python
"""
analysis of the resonance behaviour of a piezo resonator

here I will try to define a class that not only fits an admittance or impedance
circle, it shall as well be able to derive all the different frequencies
(fm, fs, fr, and fa, fp, fn) and convenience functions to calculate Qe, Qm and
all kinds of secondary characterisation data, i.e. it shall be able to auto-
matically interprete Z- and Y-circles, "the rosetta stones of transducer
analysis"

Markus Stokmaier, KIT, IKET, January 2014
"""
from os.path import join
from cPickle import Pickler
import numpy as np
from numpy import pi, sin, cos, sqrt, mean, std, amin, amax, diff
from numpy import linspace, argmin, argmax, fabs, zeros, angle, absolute, arange, array, where
from pylab import plt, cm
from scipy.signal import butter, filtfilt
from scipy.optimize import fmin, fmin_powell
from uncertainties import ufloat, umath
#import shack_postprocessing_routines as spr

def peak_sides(dat,thr=0.3,nwin=3):
    """a simple function roughly determining where the value has decayed to thr
    times the maximum value and returning the two indices"""
    mi=argmax(dat)
    mark=where(dat>thr*amax(dat),1,0)
    li,ri = None,None
    for i in xrange(mi,nwin,-1):
        if np.sum(mark[i-nwin:i])==0:
            li=i; break
    for i in xrange(mi,len(dat)-nwin):
        if np.sum(mark[i:i+nwin])==0:
            ri=i; break
    return li, ri

def is_monotonic(x):
    """found here: http://stackoverflow.com/questions/4983258/python-how-to-check-list-monotonicity"""
    dx = diff(x)
    return np.all(dx <= 0) or np.all(dx >= 0)

def rms(data,target):
    return sqrt(mean((data-target)**2))

class ZYCircle(object):

    def __init__(self,sdc,ZorY,finter=None,fakeloss=0.,afthr=0.3):
        self.sdc=sdc
        self.ZorY=ZorY
        self.finter=finter
        self.afthr=afthr # threshold used if finter='auto'
        self.fakeloss=fakeloss
        self.N=0      # number of data points
        self.f=None
        self.x=None   # x-positions of data set   --> either Z.real or Y.real
        self.y=None   # y-positions of data set   --> either Z.imag or Y.imag
        self.z=None   # x+iy
        self.Z=None
        self.R=None
        self.X=None
        self.Y=None
        self.G=None
        self.B=None
        self.phZ=None        # the angle of Z in the R-X plane
        self.phY=None        # the angle of Y in the G-B plane
        self.mgZ=None        # magnitude of Z
        self.mgY=None        # magnitude of Y
        # frequencies of impedance circle
        self.fa=None
        self.fp=None
        self.fn=None
        self.fnX=None        # where min X
        self.fmX=None        # where max X
        self.phi_a=None
        self.phi_p=None
        self.phi_n=None
        self.phi_nX=None
        self.phi_mX=None
        # frequencies of admittance circle
        self.fm=None
        self.fs=None
        self.fr=None
        self.fnB=None        # where min B
        self.fmB=None        # where max B
        self.phi_m=None
        self.phi_s=None
        self.phi_r=None
        self.phi_nB=None
        self.phi_mB=None
        self.alpha=None      # the angle with respect to the center of the fitted circle
        self.salpha=None     # slot for the smoothed angle sequence if it is decided to make it
        self.alpha_smoothed=False # whether the alpha smoothing has been applied
        self.allow_alpha_smoothing=True
        self.force_alpha_smoothing=False
        self.salpha_plotflag=True
        self.butter_order=2
        self.butter_cutoff=0.08
        self.fitcalls=0
        self.yattach=True    # whether circle has to touch imaginary axis, i.e. whether self.x == self.radius
        self.idx1=0
        self.idx2=len(self.sdc.f)
        self.ftol=1e-5
        self.ftol1d=1e-6
        self.maxiter=1000
        self.center=None     # the center of the fitted circle in the complex number plane
        self.resid=None
        self.nrms=None       # normalised root mean square distances from circle, where normalisation means division by radius
        self.offset_mean=None  # mean distance of data points from fit circle
        self.offset_std=None   # standard deviation of distance of data points from fit circle
        self.radius=None
        self.show=False      # plt.savefig or plt.show
        self.initialize()
        self.aqd={}          # analysis quantities dictionary
        self.aqdb={}         # analysis quantities dictionary
        #self.aqdc={}         # analysis quantities dictionary
        self.basic_keys=['ZorY','N','fitcalls','center','radius','resid','nrms','offset_mean','offset_std']
        self.basic_keys+=['f[0]','f[-1]','fmX','fa','fp','fn','fnX','fmB','fm','fs','fr','fnB']
        self.basic_keys+=['notes']
        self.notes=''

    def initialize(self):
        if self.finter is not None:
            if self.finter=='auto':
                if self.ZorY == 'Z':
                    redat=self.sdc.edat['Z'].real
                elif self.ZorY == 'Y':
                    redat=self.sdc.edat['Y'].real
                self.idx1,self.idx2=peak_sides(redat,self.afthr)
            else:
                self.idx1=argmin(fabs(self.sdc.f-self.finter[0]))
                self.idx2=argmin(fabs(self.sdc.f-self.finter[1]))
        self.f=self.sdc.f[self.idx1:self.idx2]
        self.N=len(self.f)
        self.Z=array(self.sdc.edat['Z'][self.idx1:self.idx2],copy=1)
        self.Y=array(self.sdc.edat['Y'][self.idx1:self.idx2],copy=1)
        self.R,self.X = self.Z.real,self.Z.imag
        self.G,self.B = self.Y.real,self.Y.imag
        self.phZ;self.mgZ = angle(self.Z),absolute(self.Z)
        self.phY;self.mgY = angle(self.Y),absolute(self.Y)
        if self.ZorY == 'Z':
            self.z,self.x,self.y = self.Z,self.R,self.X
        elif self.ZorY == 'Y':
            self.z,self.x,self.y = self.Y,self.G,self.B
        else:
            raise ValueError('wrong value for attribute self.ZorY: {}'.format(self.ZorY))

    def threshold_indices(self,thr):
        if self.ZorY == 'Z':
            mxR=amax(self.R)
            mxi=argmin(fabs(self.R-mxR))
            idx1=0; idx2=-1
            for i in xrange(mxi,-1,-1):
                if self.R[i] < thr*mxR:
                    idx1=i+1
                    break
            for i in xrange(mxi,self.N):
                if self.R[i] < thr*mxR:
                    idx2=i
                    break
        elif self.ZorY == 'Y':
            mxG=amax(self.G)
            mxi=argmin(fabs(self.G-mxG))
            idx1=0; idx2=-1
            for i in xrange(mxi,-1,-1):
                if self.G[i] < thr*mxG:
                    idx1=i+1
                    break
            for i in xrange(mxi,self.N):
                if self.G[i] < thr*mxG:
                    idx2=i
                    break
        return idx1,idx2

    def data_thinning(self,threshs,takes,throws,path,newsuffix):
        """
        example: threshs=[0.2,0.5], takes=[1,3], throws=[4,1]
        then the routine finds the stretch in the middle where G/Gmax is over the threshold, first 0.2; next, along that
        whole stretch the marker is set so that 1 data points is taken and 4 dumped; then the narrower stretch where G/Gmax
        is above 0.5 is dealt with, and there the marker is set so that 3 are taken until 1 is dumped; this of course
        overwrites part of the result of the first loop
        """
        mark=zeros(self.N,dtype=int)
        for k,thr in enumerate(threshs):
            idx1,idx2=self.threshold_indices(thr)
            take=takes[k]; throw=throws[k]; bunch=take+throw
            for i,j in enumerate(range(idx1,idx2)):
                mark[j]=int(i%bunch<take)
        self.sdc.write_selected_data(mark,path,newsuffix)

    def data_cutting(self,finter,path,newsuffix):
        self.sdc.data_cutting(finter,path,newsuffix)

    def compute_offsets(self,coords):
        self.fitcalls+=1
        xc,yc=coords
        dists=self.distances_from(xc,yc)
        if self.yattach:
            return dists - xc
        else:
            return dists - mean(dists)

    def distances_from(self,xx,yy):
        return sqrt((self.x-xx)**2 + (self.y-yy)**2)

    def residuum(self,coords):
        self.fitcalls+=1
        xc,yc=coords
        dists=self.distances_from(xc,yc)
        if self.yattach:
            return np.sum((dists - xc)**2)
        else:
            return np.sum((dists - mean(dists))**2)

    def normRMS(self):
        c,r = self.center,self.radius
        z = self.z
        return sqrt(mean((absolute(z-c) - r)**2)) / r

    def best_radius(self,point):
        xc,yc=point
        if self.yattach:
            return xc
        else:
            return mean(self.distances_from(xc,yc))

    def fit_circle(self):
        if self.ZorY == 'Z':
            center_estimate = mean(self.R), mean(self.X)
            result = fmin_powell(self.residuum, center_estimate, ftol=self.ftol, maxiter=self.maxiter, full_output=1)
            center, self.resid, direc, itr, calls, wflag =result
            self.radius=self.best_radius(center)
        elif self.ZorY == 'Y':
            center_estimate = mean(self.G), mean(self.B)
            result = fmin_powell(self.residuum, center_estimate, ftol=self.ftol, maxiter=self.maxiter, full_output=1)
            center, self.resid, direc, itr, calls, wflag =result
            self.radius=self.best_radius(center)
        self.center=center[0]+1j*center[1]
        self.nrms=self.normRMS()
        self.offset_mean = mean(fabs(absolute(self.z-self.center)-self.radius))
        self.offset_std = std(fabs(absolute(self.z-self.center)-self.radius))
        print '{} calls needed by fmin_powell to reach residuum of {} and nRMS of {}'.format(self.fitcalls,self.resid,self.nrms)
        print 'center:   x = {} and y = {}'.format(self.center.real,self.center.imag)
        print 'radius:   R = {}'.format(self.radius)
        self.update_results_dict()

    def introduce_fake_loss(self,fakeloss):
        self.fakeloss=fakeloss
        if self.ZorY == 'Z':
            self.y-=fakeloss*self.radius
        elif self.ZorY == 'Y':
            self.y+=fakeloss*self.radius

    def update_alpha(self):
        self.alpha=zeros(self.N)
        for i,z in enumerate(self.z):
            self.alpha[i]=angle(z-self.center)
        self.iron_alpha()
        if not is_monotonic(self.alpha):
            if self.allow_alpha_smoothing:
                msg='\nWarning: alpha is not monotonical!!!!!!\n'
                msg+='smooth filtered angle sequence used for determination of characteristic frequencies\n'
                msg+='better check again for validity!\n'
                self.notes+=msg
                print 40*'-'+msg+'\n'  #+40*'-'+'\n'
                # maybe the use of luxurious_smooth_finish() would be better, because it doesn't underestimate peaks
                #self.salpha=spr.luxurious_smooth_finish(self.alpha,nwin=24,offsteps=5,bds=[-4*pi,4*pi],order=3,win='hamming')
                # but because of calculation speed and the appearing smallness of the difference it makes
                # it has been decided here to go with filtfilt
                b, a = butter(self.butter_order, self.butter_cutoff)
                #self.salpha = filtfilt(b, a, self.alpha, padlen=10, padtype='odd')
                self.salpha = filtfilt(b, a, self.alpha)
                self.alpha_smoothed=True
                r=rms(self.salpha,self.alpha)
                msg='RMS difference between raw and smoothed angles: {}\n'.format(r)
                self.notes+=msg
                print msg+'\n'+40*'-'+'\n'
                if self.salpha_plotflag: self.plot_smoothed_alpha_comparison(r)
                #if r > 0.05:
                #    raise ValueError('This is an error for safety reasons: the smoothed angle data is not that close to the raw data. Better check.')
            else:
                msg='\nWarning: alpha is not monotonical!!!!!!\n'
                msg+='smoothing was not allowed\n'
                self.notes+=msg
                print 40*'-'+msg+'\n'+40*'-'+'\n'
        elif self.force_alpha_smoothing:
            b, a = butter(self.butter_order, self.butter_cutoff)
            #self.salpha = filtfilt(b, a, self.alpha, padlen=10, padtype='odd')
            self.salpha = filtfilt(b, a, self.alpha)
            self.alpha_smoothed=True
            r=rms(self.salpha,self.alpha)
            if self.salpha_plotflag: self.plot_smoothed_alpha_comparison(r)


    def iron_alpha(self):
        """roll out beyond the +pi and -pi ends, because self.w2f can't handle
        if there is a sign switch instead, because it needs an ordered self.alpha sequence"""
        idx=argmin(fabs(self.alpha))
        if is_monotonic(self.alpha):
            if self.alpha[idx+1]>self.alpha[idx]: neg=-1 # i.e. circle is in the wrong direction (it should be clockwise)
            else: neg=1
            for i in range(idx-2,-1,-1):
                if self.alpha[i]<0: self.alpha[i]+=2*pi * neg
            for i in range(idx+2,self.N):
                if self.alpha[i]>0: self.alpha[i]+=-2*pi * neg
        else:
            if mean(self.alpha[idx+1:idx+11])>mean(self.alpha[idx-10:idx]): neg=-1 # i.e. circle is in the wrong direction (it should be clockwise)
            else: neg=1
            b, a = butter(2, 0.2)
            #tmp_salpha = filtfilt(b, a, self.alpha, padlen=10, padtype='odd')
            tmp_salpha = filtfilt(b, a, self.alpha)
            mn,mx=argmax(tmp_salpha),argmin(tmp_salpha); lo,hi=min(mn,mx),max(mn,mx)
            for i in range(lo,-1,-1):
                if self.alpha[i]<0: self.alpha[i]+=2*pi * neg
            for i in range(hi,self.N):
                if self.alpha[i]>0: self.alpha[i]+=-2*pi * neg

    def w2f(self,w):
        """angle to frequency; angle is with respect to fitted circle center"""
        # the important thing: the data set usually goes around the circle in anticlockwise manner
        f=self.f
        if self.alpha_smoothed: a=self.salpha
        else: a=self.alpha
        #assert len(a)==self.N
        o1=where(a<w,1,0); s=np.sum(o1)
        #wd=360.*w/(2*pi)
        #print 'w,s: ',w,s
        assert (s!=0) and (s!=self.N)
        if o1[0]==1: j=s
        else: j=self.N-s
        i=j-1
        ratio = (w-a[i])/(a[j]-a[i])
        fval = f[i] + ratio * (f[j]-f[i])
        #print 'angle {:.1f}  {:.3f}; N {}; s {}; f[s] {}; foundf {:.1f}'.format(wd,w,self.N,s,self.f[s],fval)
        return fval

    def find_frequencies(self,return_objectives=False):
        # alternative to finding the maximum magnitude with fmin below:
        # go from 0+0j to the circle center and then extend the vector by 1*radius
        self.update_alpha()
        c=self.center; r=self.radius
        if self.ZorY == 'Z':
            obj1=lambda phi: -absolute(c+r*cos(phi)+1j*r*sin(phi))
            result1 = fmin(obj1, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_n, resid, itr, calls, wflag =result1
            obj2=lambda phi: -(c+r*cos(phi)+1j*r*sin(phi)).real
            result2 = fmin(obj2, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_p, resid, itr, calls, wflag =result2
            obj3=lambda phi: fabs((c+r*cos(phi)+1j*r*sin(phi)).imag)
            result3 = fmin(obj3, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_a, resid, itr, calls, wflag =result3
            self.phi_nX=-0.5*pi
            self.phi_mX=0.5*pi
            self.fn=self.w2f(self.phi_n[0])
            self.fp=self.w2f(self.phi_p[0])
            self.fa=self.w2f(self.phi_a[0])
            self.fnX=self.w2f(self.phi_nX)
            self.fmX=self.w2f(self.phi_mX)
            if fabs(c.imag)>r:
                self.fa=None
                self.notes+='\ncircle does not intersect with real axis\n'
        elif self.ZorY == 'Y':
            obj1=lambda phi: -absolute(c+r*cos(phi)+1j*r*sin(phi))
            result1 = fmin(obj1, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_m, resid, itr, calls, wflag =result1
            obj2=lambda phi: -(c+r*cos(phi)+1j*r*sin(phi)).real
            result2 = fmin(obj2, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_s, resid, itr, calls, wflag =result2
            obj3=lambda phi: fabs((c+r*cos(phi)+1j*r*sin(phi)).imag)
            result3 = fmin(obj3, 0., ftol=self.ftol1d, maxiter=1000, full_output=1)
            self.phi_r, resid, itr, calls, wflag =result3
            self.phi_nB=-0.5*pi
            self.phi_mB=0.5*pi
            self.fm=self.w2f(self.phi_m[0])
            self.fs=self.w2f(self.phi_s[0])
            self.fr=self.w2f(self.phi_r[0])
            self.fnB=self.w2f(self.phi_nB)
            self.fmB=self.w2f(self.phi_mB)
            if fabs(c.imag)>r:
                self.fr=None
                self.notes+='\ncircle does not intersect with real axis\n'
        self.update_results_dict()
        if return_objectives:
            if self.ZorY == 'Z':
                return obj3,obj2,obj1
            elif self.ZorY == 'Y':
                return obj1,obj2,obj3

    def pull_missing_frequencies(self,otherzyc):
        if self.ZorY == 'Z':
            self.fm,self.fs,self.fr = otherzyc.fm,otherzyc.fs,otherzyc.fr
            self.fnB,self.fmB = otherzyc.fnB,otherzyc.fmB
        elif self.ZorY == 'Y':
            self.fa,self.fp,self.fn = otherzyc.fa,otherzyc.fp,otherzyc.fn
            self.fnX,self.fmX = otherzyc.fnX,otherzyc.fmX
        self.update_results_dict()


    def to_be_annotated(self,mode):
        m=argmax(self.x)
        hi,lo = argmax(self.y),argmin(self.y)
        if mode=='triple':
            return [hi,m,lo]
        elif mode=='onlyring':
            if (m-hi > 3) and (lo-m > 3):
                return [0]+range(hi,m-1,2)+range(m-1,m+2)+range(m+3,lo+1,2)+[self.N-1]
            else:
                return [0,hi]+range(m-1,m+2)+[lo,self.N-1]
        elif mode=='quarters':
            return [0,hi,(hi+m)/2,m,(m+lo)/2,lo]

    def annotation_offsets(self,idxlist,factor=0.1,xshift=0.15):
        N=len(idxlist); offsets=zeros(N,dtype=complex)
        #print idxlist
        for i,j in enumerate(idxlist):
            #print i,j
            point=self.x[j]+1j*self.y[j]; c=self.center
            #print i,j,point,c,type(factor*(1+1.5*fabs(cos(angle(point-c))))*(point-c))
            offsets[i]=factor*(1+1.5*fabs(cos(angle(point-c))))*(point-c)
        return offsets

    def plot_overview(self,suffix=''):
        x=self.x; y=self.y; r=self.radius; cx,cy=self.center.real,self.center.imag
        ax=plt.axes()
        plt.scatter(x,y, marker='o', c='b', s=40)
        plt.axhline(y=0,color='grey', zorder=-1)
        plt.axvline(x=0,color='grey', zorder=-2)
        t=linspace(0,2*pi,201)
        circx=r*cos(t) + cx
        circy=r*sin(t) + cy
        plt.plot(circx,circy,'g-')
        plt.plot([cx],[cy],'gx',ms=12)
        if self.ZorY == 'Z':
            philist,flist=[self.phi_a,self.phi_p,self.phi_n],[self.fa,self.fp,self.fn]
        elif self.ZorY == 'Y':
            philist,flist=[self.phi_m,self.phi_s,self.phi_r],[self.fm,self.fs,self.fr]
        for p,f in zip(philist,flist):
            if f is not None:
                xpos=cx+r*cos(p); ypos=cy+r*sin(p); xos=0.2*(xpos-cx); yos=0.2*(ypos-cy)
                plt.plot([0,xpos],[0,ypos],'co-')
                ax.annotate('{:.3f} Hz'.format(f), xy=(xpos,ypos),  xycoords='data',
                            xytext=(xpos+xos,ypos+yos), textcoords='data', #textcoords='offset points',
                            arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10)
                            )
        #plt.xlim(0,0.16)
        #plt.ylim(-0.1,0.1)
        plt.axis('equal')
        if self.ZorY == 'Z':
            plt.xlabel(r'resistance $R$ in Ohm'); plt.ylabel(r'reactance $X$ in Ohm')
        if self.ZorY == 'Y':
            plt.xlabel(r'conductance $G$ in Siemens'); plt.ylabel(r'susceptance $B$ in Siemens')
        plt.title("fitting the admittance circle with Powell's method")
        tx1='best fit (fmin_powell):\n'
        tx1+='center at G+iB = {:.5f} + i*{:.8f}\n'.format(cx,cy)
        tx1+='radius = {:.5f};  '.format(r)
        tx1+='residue: {:.2e}'.format(self.resid)
        txt1=plt.text(-r,cy-1.1*r,tx1,fontsize=8,ha='left',va='top')
        txt1.set_bbox(dict(facecolor='gray', alpha=0.25))
        idxlist=self.to_be_annotated('triple')
        ofs=self.annotation_offsets(idxlist,factor=0.1,xshift=0.15)
        for i,j in enumerate(idxlist):
            xpos,ypos = x[j],y[j]; xos,yos = ofs[i].real,ofs[i].imag
            ax.annotate('{:.1f} Hz'.format(self.f[j]), xy=(xpos,ypos),  xycoords='data',
                        xytext=(xpos+xos,ypos+yos), textcoords='data', #textcoords='offset points',
                        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10)
                        )
        if self.show: plt.show()
        else: plt.savefig(join(self.sdc.plotpath,'c{}_fitted_{}_circle'.format(self.sdc.case,self.ZorY)+suffix+'.png'), dpi=240)
        plt.close()

    def plot_overview_B(self,suffix='',ansize=8,anspread=0.15,anmode='quarters',datbg=True,datbgsource=None,checkring=False):
        self.start_plot()
        if datbg: # data background desired
            self.plot_bg_data(datbgsource=datbgsource)
        #self.plot_data()
        self.plot_fitcircle()
        if checkring:
            self.plot_checkring()
        idxlist=self.to_be_annotated(anmode)
        self.annotate_data_points(idxlist,ansize,anspread)
        self.plot_characteristic_freqs(annotate=True,size=ansize,spread=anspread)
        if self.show: plt.show()
        else: plt.savefig(join(self.sdc.plotpath,'c{}_fitted_{}_circle'.format(self.sdc.case,self.ZorY)+self.sdc.suffix+self.sdc.outsuffix+suffix+'.png'), dpi=240)
        plt.close()

    def show_plot(self):
        plt.show()

    def save_plot(self,suffix=''):
        if self.center is None:
            self.fig.savefig(join(self.sdc.plotpath,'c{}_{}_circle_data'.format(self.sdc.case,self.ZorY)+self.sdc.suffix+self.sdc.outsuffix+suffix+'.png'), dpi=240)
        else:
            self.fig.savefig(join(self.sdc.plotpath,'c{}_fitted_{}_circle'.format(self.sdc.case,self.ZorY)+self.sdc.suffix+self.sdc.outsuffix+suffix+'.png'), dpi=240)
        self.fig=None; self.ax=None; plt.close()

    def clear_plot(self):
        plt.close()

    def pure_data_plot(self,connect=False,suffix='',cmap=cm.jet,bg=cm.bone(0.3)):
        #fig=plt.figure()
        ax=plt.axes()
        plt.axhline(y=0,color='grey', zorder=-1)
        plt.axvline(x=0,color='grey', zorder=-2)
        if cmap is None:
            if connect: ax.plot(self.x,self.y, 'b-',lw=2,alpha=0.5)
            ax.scatter(self.x,self.y, marker='o', c='b', s=40)
        else:
            if connect:
                if cmap in [cm.jet,cm.brg]:
                    ax.plot(self.x,self.y, 'c-',lw=2,alpha=0.5,zorder=-1)
                else:
                    ax.plot(self.x,self.y, 'b-',lw=2,alpha=0.5)
            c=[cmap((f-self.f[0])/(self.f[-1]-self.f[0])) for f in self.f]
            #c=self.f
            ax.scatter(self.x, self.y, marker='o', c=c, edgecolors=c, zorder=True, s=40) #, cmap=cmap)
        #plt.axis('equal')
        ax.set_xlim(xmin=-0.2*amax(self.x), xmax=1.2*amax(self.x))
        ax.set_aspect('equal')  #, 'datalim')
        if cmap in [cm.jet,cm.brg]:
            ax.set_axis_bgcolor(bg)
        if self.ZorY == 'Z':
            plt.xlabel(r'resistance $R$ in Ohm'); plt.ylabel(r'reactance $X$ in Ohm')
        if self.ZorY == 'Y':
            plt.xlabel(r'conductance $G$ in Siemens'); plt.ylabel(r'susceptance $B$ in Siemens')
        if self.show: plt.show()
        else: plt.savefig(join(self.sdc.plotpath,'c{}_{}_circle_data'.format(self.sdc.case,self.ZorY)+self.sdc.suffix+self.sdc.outsuffix+suffix+'.png'), dpi=240)
        plt.close()

    def start_plot(self,w=1.3,connect=False):
        self.fig=plt.figure()
        self.ax=plt.axes()
        plt.axhline(y=0,color='grey', zorder=-1)
        plt.axvline(x=0,color='grey', zorder=-2)
        self.plot_data(connect=connect)
        #plt.axis('equal')
        self.ax.set_aspect('equal', 'datalim')
        if self.center is not None:
            cx,cy=self.center.real,self.center.imag; r=self.radius
            self.ax.axis([cx-w*r,cx+w*r,cy-w*r,cy+w*r])
        else:
            xmx=amax(self.x); ymn,ymx=amin(self.y),amax(self.y)
            cx=0.5*xmx; cy=0.5*(ymn+ymx); r=0.5*(ymx-ymn)
            self.ax.axis([cx-w*r,cx+w*r,cy-w*r,cy+w*r])
        if self.ZorY == 'Z':
            plt.xlabel(r'resistance $R$ in Ohm'); plt.ylabel(r'reactance $X$ in Ohm')
        if self.ZorY == 'Y':
            plt.xlabel(r'conductance $G$ in Siemens'); plt.ylabel(r'susceptance $B$ in Siemens')

    def plot_data(self,connect=False):
        if connect: self.ax.plot(self.x,self.y, 'b-',lw=2,alpha=0.5)
        self.ax.scatter(self.x,self.y, marker='o', c='b', s=40)

    def plot_fitcircle(self):
        cx,cy=self.center.real,self.center.imag; r=self.radius
        t=linspace(0,2*pi,201)
        circx=r*cos(t) + cx
        circy=r*sin(t) + cy
        self.ax.plot(circx,circy,'g-')
        self.ax.plot([cx],[cy],'gx',ms=12)

    def plot_checkring(self):
        """
        scepticism: better check back whether the analysis results make sense, therefore let's
        proove to ourself that the equivalent circuit recreates the fitted Y- and Z-circles
        """
        w=2*pi*self.f; i=1j
        R,L,C,C0 = self.aqd['R'],self.aqd['L'],self.aqd['C'],self.aqd['C01']
        Y = i*w*C0 + 1./( R + i*w*L + 1./(i*w*C))
        Z=1./Y
        if self.ZorY=='Y':
            self.ax.plot(Y.real,Y.imag,'r--')
        if self.ZorY=='Z':
            self.ax.plot(Z.real,Z.imag,'r--')

    def annotate_data_points(self,indices,size=8,spread=0.15):
        ofs=self.annotation_offsets(indices,factor=spread,xshift=0.15)
        for i,j in enumerate(indices):
            xpos,ypos = self.x[j],self.y[j]; xos,yos = ofs[i].real,ofs[i].imag
            self.ax.annotate('{:.1f} Hz'.format(self.f[j]), xy=(xpos,ypos),  xycoords='data',
                             xytext=(xpos+xos,ypos+yos), textcoords='data', #textcoords='offset points',
                             arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10), size=size
                             )

    def plot_characteristic_freqs(self,annotate=True,size=8,spread=0.15):
        cx,cy=self.center.real,self.center.imag; r=self.radius
        if self.ZorY == 'Z':
            philist,flist=[self.phi_a,self.phi_p,self.phi_n],[self.fa,self.fp,self.fn]
        elif self.ZorY == 'Y':
            philist,flist=[self.phi_m,self.phi_s,self.phi_r],[self.fm,self.fs,self.fr]
        for p,f in zip(philist,flist):
            if f is not None:
                xpos=cx+r*cos(p); ypos=cy+r*sin(p); xos=spread*(xpos-cx); yos=spread*(ypos-cy)
                plt.plot([0,xpos],[0,ypos],'co-')
                if annotate:
                    self.ax.annotate('{:.3f} Hz'.format(f), xy=(xpos,ypos),  xycoords='data',
                                     xytext=(xpos+xos,ypos+yos), textcoords='data', #textcoords='offset points',
                                     arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10), size=size
                                     )

    def plot_bg_data(self,datbgsource=None,connect=False):
        if datbgsource is not None:
            if connect: self.ax.plot(datbgsource.edat[self.ZorY].real,datbgsource.edat[self.ZorY].imag,linestyle='-',color='grey',zorder=-100,lw=1.5,alpha=0.5)
            self.ax.scatter(datbgsource.edat[self.ZorY].real,datbgsource.edat[self.ZorY].imag,marker='x',color='grey',zorder=-100)
        else:
            if connect: self.ax.plot(self.sdc.edat[self.ZorY].real,self.sdc.edat[self.ZorY].imag,linestyle='-',color='grey',zorder=-100,lw=1.5,alpha=0.5)
            self.ax.scatter(self.sdc.edat[self.ZorY].real,self.sdc.edat[self.ZorY].imag,marker='x',color='grey',zorder=-100)

    def plot_smoothed_alpha_comparison(self,rmsval,suffix=''):
        plt.plot(self.f,self.alpha,'ko',label='data set')
        plt.plot(self.f,self.salpha,'c-',lw=2,label='smoothed angle $\phi$')
        plt.xlabel('frequency in Hz')
        plt.ylabel('angle $\phi$ in coordinates of circle')
        plt.legend()
        ylims=plt.axes().get_ylim()
        plt.yticks((arange(9)-4)*0.5*pi, ['$-2\pi$','$-3\pi/2$','$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'])
        plt.ylim(ylims)
        plt.title('RMS offset from smooth curve: {:.4f}'.format(rmsval))
        if self.show: plt.show()
        else: plt.savefig(join(self.sdc.plotpath,'salpha','c{}_salpha_on_{}_circle'.format(self.sdc.case,self.ZorY)+self.sdc.suffix+self.sdc.outsuffix+suffix+'.png'), dpi=240)
        plt.close()

    def analyse_circle(self):
        if self.ZorY=='Z':
            self.analyse_impedance_circle()
        elif self.ZorY=='Y':
            self.analyse_admittance_circle()

    def analyse_impedance_circle(self):
        # all quantities for equivalent circuit components have a name like self.ec_*
        d=self.aqd; db=self.aqdb; #dc=self.aqdc
        d['Rmax']=2*self.center.real
        d['X0']=self.center.imag
        if (self.fr is not None) and (self.fa is not None):
            db['Qprod'] = self.fr**2 / (self.fa**2 - self.fr**2)
        if (self.fs is not None) and (self.fp is not None):
            db['k'] = umath.sqrt(1 - self.fs**2/self.fp**2)
        d['header']="1st run based on own target resonance"
        db['header']="2nd run incorporating fr and fs"
        #dc['header']=""

    def analyse_admittance_circle(self):
        # all quantities for equivalent circuit components have a name like self.ec_*
        d=self.aqd; db=self.aqdb; #dc=self.aqdc
        d['Gmax']=2*self.center.real
        d['R']=1./d['Gmax']
        d['Bs']=self.center.imag
        d['Qe']=d['Bs']/d['Gmax']
        d['Qm'] = self.fs / (self.fnB-self.fmB)
        d['Qprod']=d['Qe']*d['Qm']
        d['L'] = d['Qm'] * d['R'] / (2*pi*self.fs)
        d['C'] = 1. / (d['Qm']*d['R']*2*pi*self.fs)
        d['C01'] = d['Bs'] / (2*pi*self.fs)
        d['C02'] = d['C'] * d['Qprod']    # based on Bs, Gmax, fs, fnB, fmB
        d['k'] = umath.sqrt(d['C']/(d['C01']+d['C']))
        d['Ma'] = d['Qm'] / d['Qprod']
        d['Mb'] = 1. / (2*pi*self.fs*d['C01']*d['R'])
        d['Gamma'] = d['C']*3e-3 / (2*pi*34e-3*25e-3)
        if (self.fs is not None) and (self.fp is not None):
            db['k'] = umath.sqrt(1 - self.fs**2/self.fp**2)
            db['Qprod2'] = (1-db['k']**2) / db['k']**2
            db['M2'] = d['Qm'] / db['Qprod2']
        if (self.fa is not None) and (self.fr is not None):
            db['Qprod1'] = self.fr**2 / (self.fa**2 - self.fr**2)
            db['M1'] = d['Qm'] / db['Qprod1']
            db['C01'] = d['C'] * db['Qprod1'] # based on fa, fr
            db['C02'] = d['C'] * db['Qprod2'] # based on fs, fp
            db['k1'] = umath.sqrt(d['C']/(d['C01']+d['C']))
            db['k2'] = umath.sqrt(d['C']/(d['C02']+d['C']))
            db['k1b'] = umath.sqrt(d['C']/(db['C01']+d['C']))
            db['k2b'] = umath.sqrt(d['C']/(db['C01']+d['C']))
            #db['Bs'] = 2*pi*self.fs * db['C0']
            #db['Qe'] = db['Qprod'] / d['Qm']
        #d['Qprodb'] = d['C0']/d['C']
        d['header']="1st run based on own target resonance"
        db['header']="2nd run incorporating fa and fp"
        #dc['header']="3rd run with alternative formulas"

    def update_results_dict(self):
        for key in self.basic_keys:
            self.aqd[key]=eval('self.'+key)

    def analysis_printout(self):
        ofl=open(join(self.sdc.textpath,self.ZorY+'_circle_ana_c'+str(self.sdc.case).zfill(2)+'_'+self.ZorY+self.sdc.suffix+self.sdc.outsuffix+'.txt'),'w')
        for word in ['sdc.case','sdc.suffix','sdc.datapath']:
            txt='{} :  {}\n'.format(word,eval('self.'+word))
            ofl.write(txt); #print txt
        ofl.write('\n\n')
        self.write_down_dict(self.aqd,ofl)
        self.write_down_dict(self.aqdb,ofl)
        #self.write_down_dict(self.aqdc,ofl)
        ofl.write(self.notes)
        ofl.close()

    def write_down_dict(self,ddict,ffile):
        ffile.write(ddict['header']+'\n')
        keylist=ddict.keys()
        keylist.sort()
        for word in keylist:
            if word not in ['notes','header']:
                txt='{} :  {}\n'.format(word,ddict[word])
                ffile.write(txt); #print txt
        ffile.write('\n\n')


    def introduce_uncertainty(self,attr,u):
        #uf=ufloat(str(eval('self.'+attr))+'('+str(u)+')')
        #exec('self.'+attr+'=ufloat('+str(eval('self.'+attr))+'('+str(u)+'))')
        val=eval('self.'+attr)
        exec('self.'+attr+'=ufloat("{}({})")'.format(val,u))

    def pickle_aqd(self,suffix=''):
        ofile=open(join(self.sdc.picklepath,'zyc_aqd_pickle_c'+str(self.sdc.case).zfill(2)+'_'+self.ZorY+self.sdc.suffix+self.sdc.outsuffix+suffix+'.txt'), 'w')
        einmachglas=Pickler(ofile)
        einmachglas.dump(self.aqd)
        ofile.close()

    def pickle_self(self,suffix=''):
        ofile=open(join(self.sdc.picklepath,'zyc_pickle_c'+str(self.sdc.case).zfill(2)+'_'+self.ZorY+self.sdc.suffix+self.sdc.outsuffix+suffix+'.txt'), 'w')
        einmachglas=Pickler(ofile)
        einmachglas.dump(self)
        ofile.close()
