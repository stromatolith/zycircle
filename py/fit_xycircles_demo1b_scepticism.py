#!python
"""
ZYCircle demo 1

a) data is being read in the constructor of the DataContainer
b) fitting admittance circle
c) fitting impedance circle

d) scepticism: better check back whether the analysis results make sense,
   therefore let's proove to ourselves that the equivalent circuit with the
   quantities gained from analysing the admittance circle really recreates
   the fitted admittance circle (red circle in the plot)

Markus Stokmaier, KIT, IKET, March 2014
"""
import helpers as hlp
from ZYCircle import ZYCircle

dc=hlp.DataContainer(129,switchiphase=True)

yc=ZYCircle(dc,'Y',finter=[18400,18800])
yc.fit_circle()
yc.find_frequencies()
yc.analyse_circle()
yc.plot_overview_B(datbg=True,checkring=True)
yc.analysis_printout()

zc=ZYCircle(dc,'Z',finter=[18700,19200])
zc.fit_circle()
zc.find_frequencies()
zc.analyse_circle()

# now transferring stuff between the analysis quantities dictionaries
for key in ['R','L','C','C01']:
    zc.aqd[key]=yc.aqd[key]
zc.plot_overview_B(datbg=True,checkring=True)
zc.analysis_printout()



