zycircle
========

python scripts for analysing admittance and impedance circles

These are postprcessing scripts I wrote for postprocessing resonator characterisation data. I worked on a resonator driven by a piezoelectric transducer and I needed to analyse frequency sweep data long after recording far away from the lab and any network analyser or LCR meter.

####contents:
- `helpers.py` contains a class for reading data and some plotting, but also some functions for filtering and Q-factor determination
- `ZYCircle.py` contains the class for analysing admittance and impedance circles
- `simple_plot_calls.py`: demo for plotting utilities of the `DataContainer` class
- `fit_xycircles_demo1.py`: simplest use case of analysing admittance and impedance circles
- `fit_xycircles_demo2.py`: automatic frequency interval detection of the fit-relevant data subset
- `fit_xycircles_demo3.py`: data exchange between Y-circle instance and Z-circle instance, more analysis quantities can then be computed
- `luxury_smooth_demo.py`: this is a sort of filter that might be of some interest. It shows the performance of a data smoothing function I programmed (not nicely though) which has the ability of cleaning a lot of noise while not suffering from two sometimes important drawbacks: it doesn't broaden peaks like a simple lowpass, and it doesn't turn poles into wavelets.


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