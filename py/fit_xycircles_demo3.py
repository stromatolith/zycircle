#!python
"""
ZYCircle demo 1

Here comes the problem: there are important formulae where analysis quantities
from both, the Y- and the Z-circle are needed, e.g.

           k  =  1 - fs**2/fp**2

that's of course no problem here:
a) circle analysis once
b) data exchange between Z- and Y-circle
c) another round of analysis

Markus Stokmaier, KIT, IKET, March 2014
"""
import helpers as hlp
from ZYCircle import ZYCircle

dc=hlp.DataContainer(129,switchiphase=True)

yc=ZYCircle(dc,'Y',finter=[18400,18800])
yc.fit_circle()
yc.find_frequencies()
yc.plot_overview_B(datbg=True)
yc.analyse_circle()
#yc.analysis_printout()

zc=ZYCircle(dc,'Z',finter=[18700,19200])
zc.fit_circle()
zc.find_frequencies()
zc.plot_overview_B(datbg=True)
zc.analyse_circle()
#zc.analysis_printout()

yc.pull_missing_frequencies(zc) # data exchange happens here
zc.pull_missing_frequencies(yc) # data exchange happens here

yc.analyse_circle()
yc.analysis_printout()

zc.analyse_circle()
zc.analysis_printout()

"""
the only thing you have to do is to create a few folders:
    
+-workdir
  |
  +-analysis    -> the text output will land here
  |
  +-data        -> the provided test data has to be here
  |
  +-plots       -> plots are saved in here
  | |
  | +-salpha    -> an important type of check plot is saved here (salpha stands
  |                for "smoothed alpha", whereby alpha corresponds to the angle
  |                phi in the explaining text below)
  |
  +-py          -> the location of the python scripts
"""

"""
one more thing:
You may ask: why did you make the ZYCircle class so dumb? Why can it be only
one circle and not both? Why was it not made to be powerful enough to treat the
whole dataset and analyse both circles together, then there would be no need
here in the main script to manually command the information exchange?

My answer: if you did that, you would run the danger of swamping the ZYCircle
class definition with loads of code dealing with robustness side issues to make
the automatic detection, separation, and analysis of the two circles possible
without being too often mislead on poor-quality data. If you outline the task
too large, it forces you to think of so many eventualities, to anticipate so
many potential bugs, that it's too large drag. Therefore I decided to keep it
simple (or at least as simple as possible, as the ZYCircle source code already
suffers from being cluttered with the plotting convenience stuff.
"""