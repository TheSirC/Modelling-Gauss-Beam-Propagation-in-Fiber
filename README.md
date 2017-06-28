# Modelling focalisation of a gaussian beam through a cyndrical interface
![Animation](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/blob/master/Simulations/Results/animation.gif)

## Theory
The theory used to model the behavior of the beam is describded in [P. A. Bélanger, "Beam propagation and the ABCD ray matrices," Opt. Lett. 16, 196-198 (1991)](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/blob/master/Recherches/B%C3%A9langer%20-%201991%20-%20Beam%20propagation%20and%20the%20ABCD%20ray%20matrices.pdf).

## Programming 
### Intuition of the program
The code computes each plane starting from the far left to the far right (the direction of propagation through the cylinder) and computes both direction independently (in [`calcul_x`](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/blob/master/Simulations/calcul_x.m) and [`calcul_y`](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/blob/master/Simulations/calcul_y.m) respectively). For each plane, it computes the field then the intensity. After getting the two intensity planes, a threshold (found experimentally for the sillica, 1/e^2 for the first version) is given to be used as the radii (semi-minor and semi-major axis) of the ellipse (in [`calcul_ellipse`](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/blob/master/Simulations/calcul_ellipse.m)). Each ellipse is drawn to the 3D plot one after the other (again following the propagation direction).

## Details
For further details on the implementation, please refer to the [code](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/tree/master/Simulations) which is commented extensively.

## Results
You can find images (intensity and isometric view) from different angles (XZ plane, YZ planes) for different depths [here](https://github.com/TheSirC/Modelling-Gauss-Beam-Propagation-in-Fiber/tree/master/Simulations/Results/).
