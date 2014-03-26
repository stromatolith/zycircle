#!python
"""
as a demo of the DataContainer class and other functions in helpers.py a script
calling the plotting functions I used most often.

of general interest:
    - the Q-factor function: how it works with noise, telling error messages when it can't work
    - the smoothing function for phase signals (or other signals cycling through bounded domains)

The rest of the stuff is maybe not that much of general interest and perhaps
only an example of matplotlib usage.

Markus Stokmaier, KIT, IKET, March 2014
"""
from os import getcwd
from os.path import join, dirname
from numpy import array, amax
from numpy.random import randn
import matplotlib.pyplot as plt
import helpers as hlp

dca=hlp.DataContainer(95,switchiphase=True)
dcb=hlp.DataContainer(129,switchiphase=True)
dcc=hlp.DataContainer(125,switchiphase=True)

#for dc in [dca,dcb]:
#    dc.electric_rawdat_plot()
#    dc.outsuffix='_log'
#    dc.electric_rawdat_plot(logI=True,logU=True)
#    dc.outsuffix=''
#    dc.electric_plot()
#    dc.standard_plot()
#
#dcc.compare_with(dcb)
#dcb.compare_with(dca)
#dcb.compare_electric_rawdat_plot(dca,logI=True,logU=True)
#dcc.compare_electric_rawdat_plot(dcb,logI=True,logU=True)
#
#"""
#Q-factor computation
#arg 1: data type
#arg 2: frequency interval
#arg 3: window size for peak detection
#arg 4: window size for filtering
#"""
#for dc in [dca,dcb]:
#    dc.determine_Q('hyd',[19700,20100],20)
#    dc.determine_Q('mic',[19700,20100],20)
#    dc.determine_Q('hgain',[19700,20100],20)
#    dc.determine_Q('mgain',[19700,20100],20)
#dca.determine_Q('mgain',[18500,18750],20)
dca.outsuffix='_multipeaks_1'
dca.determine_Q('mgain',[dca.f[0],dca.f[-1]],20,nwin_s=8,smoothing='standard')
dca.outsuffix='_multipeaks_2'
dca.determine_Q('mgain',[dca.f[0],dca.f[-1]],20,nwin_s=8,smoothing='luxury')
#"""
#how the Q-factor detection deals with noise: Q-factor of filtered blue curve
#"""
#dc_bad=hlp.DataContainer(129,switchiphase=True)
#mic=array(dc_bad.data['mic'],copy=1)
#n=len(dc_bad.data['mic'])
#mx=amax(dc_bad.data['mic'])
#for i,level in enumerate([0.05,0.1,0.2,0.3]):
#    dc_bad.data['mic'] = mic + level*mx*randn(n)
#    dc_bad.outsuffix='_with_noise_{}'.format(i)
#    dc_bad.determine_Q('mic',[19700,20100],20)
#
#"""
#the smoothing function hlp.luxurious_smooth_finish()
#Everybody knows what happens when you smooth data having a pole, i.e. data that
#leaves the domain at the top and then comes in from below: you turn the pole
#into a wavelet! But this smoothing function doesn't do that. Also, it doesn't
#wash out small structures too much depending on the settings)
#"""
#dcd=hlp.DataContainer(296,switchiphase=True)
#freq=dcd.f
#mic_phs=dcd.data['mic_phs']
#mic_sphs=hlp.luxurious_smooth_finish(mic_phs,nwin=24,offsteps=5,bds=[-180.,180.],order=4,win='hamming')
#plt.plot(freq,mic_phs,'ko')
#plt.plot(freq,mic_sphs,'co',markersize=2.5,mec='c')
#plt.xlabel('frequency in Hertz')
#plt.ylabel('microphone phase')
#plt.yticks([-180,-135,-90,-45,0,45,90,135,180])
#plt.title('an expensive but performant smoothing function\nfor signals cycling through an interval e.g. phase signals')
#plt.savefig(join(dcd.plotpath,'luxury_smoothed_phase_signal.png'))
##plt.show()
#plt.close()
