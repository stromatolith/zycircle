#!python
"""
ZYCircle demo 1

Here comes the problem: there are important formuae where analysis quantities
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
