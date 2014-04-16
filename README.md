zycircle
========

**python scripts for analysing admittance and impedance circles**

These are postprcessing scripts I wrote for postprocessing characterisation data recorded from a resonator driven by a piezoelectric transducer. The coordinates and characteristic frequencies of admittance and impedance circles can tell you a lot about a transducer. Piezoelectric transducers can often be represented by such an equivalent circuit:
~~~~~~ python
             |
             |
       o-----o-----o
       |           |
       |           Z   R
       |           Z
       |           |
 C_0  ===          S   L
       |           S
       |           |
       |          ===  C
       |           |
       o-----o-----o
             |
             |
~~~~~
Here, the R-L-C branch represents the spring-mass system of the piezoelectric transducer in which energy sloshes back and forth between the two forms of elastic deformation and kinetic energy. The other branch with the parallel capacitance C_0 (also called parasitic capacitance) stands for the AC leak path existing because cables and electrodes act as a capacitor of some size (as I understand it so far).

Admittance circles are the "Rosetta stone of transducer analysis", because the equivalent circuit quantities can be easily determined from them. But why should one be interested in manually examining admittance and impedance circles? Why not simply use an LCR-meter or network analyser? Well, maybe
- you want to be in control of which data points are used for curve or circle fitting,
- using a network analyser you only get the end result in terms of sizes of L, C, R, and some errors or deviations but you never get to see the fitted curves or circles plotted on the actual data, and if you have noisy or otherwise poor data that visualisation is much more helpful than some error bounds for judging whether the analysis results make some sense,
- you want to backcheck the performance of your network analyser on poor data and gain more trust in it,
- you want to compare the results of using alternative formulae, or
- (my motivation is somewhat more embarassing) you learn that a look at admittance circles would make sense only long after the measurement campaign.

Where does the actual formula crunching happen?
- The class `ZYCircle` contains a method `find_frequencies()`, here the characteristic frequencies (where max real value, max magnitude, cutting real axis, etc.) are determined.
- The class `ZYCircle` has the methods `analyse_admittance_circle()` and `analyse_impedance_circle()`, here happens all the important maths. The result values are stored in python dictionaries allowing a slim printout function.


####contents:
- `helpers.py` contains a class for reading data, called `DataContainer`, and some plotting, but also some functions for filtering and Q-factor determination
- `ZYCircle.py` contains the class for analysing admittance and impedance circles
- `simple_plot_calls.py`: demo for plotting utilities of the `DataContainer` class
- `fit_xycircles_demo1.py`: simplest use case of analysing admittance and impedance circles
- `fit_xycircles_demo1b_scepticism.py`: do the gained values for R, L, C, C0 recreate the circle?
- `fit_xycircles_demo2.py`: automatic frequency interval detection of the fit-relevant data subset
- `fit_xycircles_demo3.py`: data exchange between Y-circle instance and Z-circle instance, more analysis quantities can then be computed
- `luxury_smooth_demo.py`: this is a sort of filter that might be of some general interest also outside the context of transducer analysis. The demo shows the performance of that data smoothing function I programmed (not nicely though) which has the ability of cleaning a lot of noise while not suffering from two sometimes important drawbacks: it doesn't broaden peaks like a simple lowpass, and it doesn't turn poles into wavelets.
- there are some test data sets located in the folder `data` which are needed by the demo scripts


####References:

I compiled the formulae after having found them here:

@article{ieee_std177_1966,
	title = {{IEEE} Standard Definitions and Methods of Measurement for Piezoelectric Vibrators},
	doi = {10.1109/IEEESTD.1966.120168},
	journal = {{IEEE} Std No.177},
	year = {1966},
	pages = {1--19}
},

@article{ieee_std176_1988,
	title = {{IEEE} Standard on Piezoelectricity},
	doi = {10.1109/IEEESTD.1988.79638},
	journal = {{ANSI/IEEE} Std 176-1987},
	year = {1988},
	pages = {0\_1--}
},

@book{wilson_introduction_1988,
	address = {Los Altos, {CA}},
	title = {Introduction to theory and design of sonar transducers},
	isbn = {9780932146229},
	language = {en},
	publisher = {Peninsula},
	author = {Wilson, Oscar Bryan},
	year = {1988}
},

@inproceedings{deangelis_optimizing_2010,
	title = {Optimizing piezoelectric ceramic thickness in ultrasonic transducers},
	doi = {10.1109/UIA.2010.5506062},
	booktitle = {Ultrasonic Industry Association Symposium ({UIA)}, 2010 39th Annual},
	author = {{DeAngelis}, {D.A.} and Schulze, {G.W.}},
	year = {2010},
	pages = {1--9}
}