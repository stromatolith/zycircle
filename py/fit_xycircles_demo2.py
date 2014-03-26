#!python
"""
ZYCircle demo 1

a) cutting one resonance out and saving the cloipped data
b) use finter='auto' with a threshold of 0.08
   this means: for circle fitting you don't have to think about which
   frequency interval you should best specify manually; instead, only the
   points in the data set lying above a threshold of 0.08*Gmax (or 0.08*Rmax)
   will be used, or more precisely, a corresponding frequency interval will be
   figured out (note that thus the smaller circle from the other resonance will
   be ignored as should be the case; the used points are the blue ones)
c) fitting impedance circle

Markus Stokmaier, KIT, IKET, March 2014
"""
import helpers as hlp
from ZYCircle import ZYCircle

dc=hlp.DataContainer(129,switchiphase=True)

yc=ZYCircle(dc,'Y',finter='auto',afthr=0.08)  # "afthr" stands for auto finter threshold
yc.fit_circle()
yc.find_frequencies()
yc.plot_overview_B(datbg=True)
yc.analyse_circle()
yc.analysis_printout()

zc=ZYCircle(dc,'Z',finter='auto',afthr=0.08)  # "afthr" stands for auto finter threshold
zc.fit_circle()
zc.find_frequencies()
zc.plot_overview_B(datbg=True)
zc.analyse_circle()
zc.analysis_printout()


"""
But how about the smaller resonance?
--> the frequency interval has to be specified only roughly to cut out the
right piece, then finter='auto' as above can be used again
"""
from os import getcwd, mkdir
from os.path import join, dirname

loc=dirname(getcwd())
datapath=join(loc,'custom_data')
try: mkdir(datapath)
except: pass
dc.data_cutting([19600,20100],path=datapath,newsuffix='_res2')

dc2=hlp.DataContainer(129,switchiphase=True,datapath=datapath,suffix='_res2')

yc2=ZYCircle(dc2,'Y',finter='auto',afthr=0.08)  # "afthr" stands for auto finter threshold
yc2.fit_circle()
yc2.find_frequencies()
yc2.plot_overview_B(datbg=True)
yc2.plot_overview_B(datbg=True,datbgsource=dc,suffix='_allbg')
yc2.analyse_circle()
yc2.analysis_printout()

zc2=ZYCircle(dc2,'Z',finter='auto',afthr=0.08)  # "afthr" stands for auto finter threshold
zc2.fit_circle()
zc2.find_frequencies()
zc2.plot_overview_B(datbg=True)
zc2.plot_overview_B(datbg=True,datbgsource=dc,suffix='_allbg')
zc2.analyse_circle()
zc2.analysis_printout()


