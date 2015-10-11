# Running Coupling in the Standard Model and Beyond

This directory lists the code I used as part of my investigation into running coupling.

## A Note

I am now a professional software developer. This was written in 2011, when I had no formal software training. Bear that in mind if you're looking through this for an interview or something! (As a footnote, I've seen worse code in academic circles...)

Each of the files has a header explaining roughly what is contained within. Also, most methods have comments to guide you.

## The Code

The following libraries are used:

- [Models](models.py) contains the theoretical models used in this investigation:
  - The Standard Model (`StandardModel`)
  - The Minimally Supersymmetric Standard Model (`MSSMModel`)
  - The Minimally Supersymmetric Standard Model with Extra Compactified Dimensions (`MSSMExtraDim`)
- [Methods](methods.py) contains the mathematical methods used to find the minimum value of a function in a given domain:
  - Monte Carlo minimisation - evaluate a function at random points in a multi-dimensional domain for a fixed number of trials.
  - Systematic minimisation - evaluate a function at all points in a multi-dimensional array.
- [Plots](plots.py) contains the plots output following unification:
  - Coupling Plot - shows the three fundamental forces' strengths as a function of the energy level; error bars are also plotted.
  - Parameter Contour - 2D contour plot showing how the degree of unification varies across two parameters.
  - Parameter 3D - 3D version of the Parameter Contour.

The remaining files use the libraries to perform analysis:

- [SM-MSSM](sm-mssm.py) analyses the SM and MSSM, generating Coupling Plots.
- [MSSM-BP](mssm-bp.py) implements the Monte-Carlo Levenberg-Marquardt (MCLM) algorithm for minimisation and uses it to produce contour plots. It does this using back-propagation: setting the unification point and seeing how close it gets to experimental values when returning to observed energy levels.
- [ED-MSSM](ed-mssm.py) provided a number of analyses when extra dimensions are added to the MSSM.

These files contain large amounts of commented-out code. Sorry.

## The Results

These scripts take a long time to run - hey, we're solving simultaneous differential equations here!

Check out the figures in the reports at the root of this directory to see what comes out of these scripts.
