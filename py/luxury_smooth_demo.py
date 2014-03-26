#!python
"""
I have to admit, I afforded myself a rather expensive data filter, i.e. it is
really unbelievably slow, but there are cases where you care less about speed
than about a few relevant data sets being treated the right way

Everybody knows what happens when you smooth data having a pole, i.e. data that
leaves the domain at the top and then comes in from below: you turn the pole
into a wavelet! But this smoothing function doesn't do that.

Secondly, it doesn't wash out small structures too much. Why is this helpful
for Q-factor determination? If the data is noisy you can't just go down the
sides of a peak and measure how far it goes until the first data point lies
below the 1/sqrt(2) threshold. But if you filter with a simple lowpass until
there are no more local minima along the peak sides, then you wash out the
whole peak form and destroy what you are trying to quantify, the peak shape. A
solution would be to fit some sort of bell function. The drawbacks: your
examined peak is not necessarily of the shape of your model function, and
secondly, structures at the foot of the peak or overlapping other peaks mess
up the fitting procedure even though you ideally would like to not care about
the foot or neighbourhood. Below you see my homebrew solution to the dilemma.

Markus Stokmaier, KIT, IKET, March 2014
"""
from os.path import join
import matplotlib.pyplot as plt
import helpers as hlp

#--- case 1: phase signal smoothing
dca=hlp.DataContainer(296,switchiphase=True)
freq=dca.f
mic_phs=dca.data['mic_phs']
mic_sphs=hlp.luxurious_smooth_finish(mic_phs,nwin=24,offsteps=5,bds=[-180.,180.],order=4,win='hamming')
plt.plot(freq,mic_phs,'ko')
plt.plot(freq,mic_sphs,'co',markersize=2.5,mec='c')
plt.xlabel('frequency in Hertz')
plt.ylabel('microphone phase')
plt.yticks([-180,-135,-90,-45,0,45,90,135,180])
plt.title('an expensive but performant smoothing function\nfor signals cycling through an interval e.g. phase signals')
plt.savefig(join(dca.plotpath,'luxury_smoothed_phase_signal.png'))
#plt.show()
plt.close()

#--- case 2: use case Q-factor determination
dcb=hlp.DataContainer(95,switchiphase=True)
dcb.outsuffix='_multipeaks_1'
dcb.determine_Q('mgain',[dcb.f[0],dcb.f[-1]],20,nwin_s=8,smoothing='standard')
dcb.outsuffix='_multipeaks_2'
dcb.determine_Q('mgain',[dcb.f[0],dcb.f[-1]],20,nwin_s=8,smoothing='luxury')


"""
how does the function work?

a) take a piece of the data[i:i+n]
b) plynomial fitting (special treatment: rolling data upwards or downwards if
   it leaves and reenters the domain on the other side)
c) i+=1

d) cycle a-b-c until you're through the data

e) optionally multiply each piece with a window, e.g. Hamming
f) at the end put all the pieces together again
"""
