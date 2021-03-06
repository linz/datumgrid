datumgrid
=========

Program to build a grid defining a one or two dimensional datum offset based on values observed
at a set of control points.  

This uses an approach of directly calculating the grid offset using a least squares calculation 
rather than interpolating control point observations onto the grid nodes.  This approach is preferred
where the datum offset reflects historical systematic errors such as in levelling datums or historical
horizontal datums.  There is no good basis for statistical interpolation such as inverse distance 
weighted schemes.  On the other hand there are some expectations on the generated grid, such
as not introducing unnecessary distortion between the datums.

The approach used is to solve using least squares and including "observations" to minimize two quantities,
firstly the fit of the observations to the control points, and secondly the distortion within each
grid square.  Least squares is used here not as a statistical technique but instead as a mathematically
simple optimisation technique.  Where warranted this programme can be run iteratively (eg by building 
a script using it) to emulate more robust weighting schemes.

The main disadvantage of this approach of IDW is that it does not scale well - as the 
grid size increases the least squares matrices become very big and solving can take a long time!
Similar results may be obtained more efficiently by kriging, especially for one dimensional data
(height offsets).  For two dimensional data the kriging method would need to be built based
on a model of affine transformations. 

The program uses two input files, one is a control point file, and the other a configuration file
defining parameters of the fit.  See the [datumgrid.help](https://raw.githubusercontent.com/linz/datumgrid/master/src/datumgrid.help) file for more information about these.

The analyse_data python script is an example of an iterative reweighting script.  The other scripts
are best ignored!

This is definitely a "research" program .. though it is only slightly modified from the version that 
was used to create the official transformation from the NZGD1949 datum to the NZGD2000 datum in 1999.
