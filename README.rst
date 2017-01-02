LASDiA: Liquids and Amorphous Solids Diffaction data Analyzer
=============================================================

LASDiA is a software developed to analyze diffraction data from liquids and
amophous solids.

LASDiA can be used installing Python and running the test (LASDiAScript.py)
or GUI (LASDiA.py) script you can find in its folder.

Dependencies
------------

LASDiA is developed in Python 3.5 in Windows and MacOs environment.
The following modules are necessary:

	* numpy == 1.11.2
	* scipy == 0.17.1
	* matplotlib == 1.5.3
	* PyQt5 == 5.7

Windows distributions
---------------------

In Windows environment, LASDiA is developed using WinPython 3.5.2.2Qt5.
This version present the following modules:

	* numpy == 1.11.1
	* matplotlib == 1.5.2

this create a compatibility problem with the embedded matplotlib widget in PyQt5.
To resolve this problem follow these steps in this exact order because matplotlib
will install a numpy version incompatible with WinPython:

	* update matplotlib using pip 
	* install the numpy mkl version, it can be find on
	  Christoph Gohlke's website: http://www.lfd.uci.edu/~gohlke/pythonlibs/

Tests with WinPython 3.5.2.3Qt5 will be done asap.

MacOs distributions
-------------------

In MacOs environment, LASDiA can be used installing the official Python 3.5 and
updating the modules with pip.