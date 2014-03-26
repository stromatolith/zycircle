#!python
"""
ZYCircle demo 1

a) data is being read in the constructor of the DataContainer
b) fitting admittance circle
c) fitting impedance circle

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
yc.analysis_printout()

zc=ZYCircle(dc,'Z',finter=[18700,19200])
zc.fit_circle()
zc.find_frequencies()
zc.plot_overview_B(datbg=True)
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


about the checkplots in the smoothed_phi folder:

Once there is a fitted Y- or X-circle, then the points of the data set can be
expressed in terms of distance and angle with respect to the circle center, i.e.
in the cylindrical coordinate system coincident with the circle. If the data
is nice, then the angle, let's call it phi, reduces from step to step. But if
the data is not so nice, if it is a bit more scattered around the circle line,
and if the current-vs-voltage phase signal had some noise in the measurement,
then the phi does not decrease monotonically.
This means that going through the data set point by point (\ie in clockwise
direction through the circle), there are intermittent forward and backward
steps. How should the characteristic frequencies ($f_r$,$f_a$,$f_{mB}$
$f_{nX}$,...) be determined in such a case? To take just the data point nearest
to the corresponding point of the fitted circle would be suboptimal. And what
should be taken as nearest, smallest distance or smallest angular offset? The
approach to choose the data point with the largest $G$ to read out $f_s$, the
one with the largest $|Z|$ for $f_n$, the smallest $X$ for $f_a$ \etc would be
even worse.

The approach followed here was to write the data set in a new coordinate system,
a cylindrical one corresponding directly to the fitted $Y$- or $Z$-circle. In
the cylindrical system each point is defined by an angle $\varphi$ and a radius
$r$, so in the case of an admittance circle one has
\begin{eqnarray*}
 G &=& C_G +r\cos \varphi \\
 B &=& C_B +r\sin \varphi
\end{eqnarray*}
where $C=(C_G,C_B)$ is the center of the fitted $Y$-circle. The checkplots
created in the smoothed_phi folder show the angle $\varphi$ over the frequency
$f$ for the data set underlying the $Z$- or $Y$-circle. Smoothing the data set
can yield a function (the cyan line) which can be used for translating $f$ into
$\varphi$ and vice versa. In practice this is done by linear interpolation of
the smoothed data set. Calculating the characteristic frequencies that way
yields much more telling values, because they are not influenced by the random
offset of single data points. For example, $f_s$ can be computed by requesting
the frequency value corresponding to $\varphi=0$ from the interpolation
function, or $f_r$ upon first calculating the angle $\varphi$ where the fitted
circle crosses the abscissa in the $G$-$B$-coordinate system and then
requesting the frequency for that angle. The lowpass filter used for the
smoothing was a Butterworth filter of order two, and the cutoff has been chosen
to be \num{0.08} times the Nyquist frequency. The filter has been applied
bidirectionally. This has been done using the functions \verb|butter| and
\verb|filtfilt| of the \verb|scipy.signal| section of the SciPy library
\cite{scipy}.

"""


