Mathematica
===========

A collection of Mathematica homebrew packages, mostly related to
Control, Dynamics, and Diff. Geom.

## Description

The set of libraries is borken down into the following types and
packages:

Name | Provides
-----| --------
Packages  | Package handling, for specifying and loading other packages
Basic   | Some basic extensions or definitions to Mathematica
Math      | Basic additions
Dynamics   | Dynamical Systems
Averaging | for time-periodic systems and averaging
DiffGeometry | Differential Geometry
Graphics   | Graphics and animation
Mechanics  | Lagrangian (and maybe Hamiltonian) Mechanics
GeoMechanics | Geometric Mechanics
Robotics    | Robotics related material

A few observations.

**Object-Oriented.**
Back when this library was implemented (2001-2003), Mathematica had
little to no support for object oriented programming.  A quick google
search seems to indicate that this still holds.  There are more options
now than before, but nothing seems to be definitive.
I kludged it in the older Mathematica versions a decade ago, so it   This
needs revision.

**Need Updating.**
They were written a decade ago and may now be superceded by actual
Mathematica packages.  I imagine, we'll figure this out as they get
used and updated.  A lot of the packages are really for manipulating
things symbolically for use elsewhere.  Nevertheless, there is also some
code to actually process or simulate things numerically within
Mathematica.  

**Numerical Integration.**
In general, I found the language design to be a bit opposed to numerical
implementations, so there was tendency to work things out enough in
Mathematica then move to Matlab.


## The packages.

### Basic

Just adds some basic definitions or functions that allows me to
customize the Mathematica experience towards the way I conceive of
things, or the way I like to see the output.



### Math

Some additions that make existing Mathematica code conform a bit more to
the written math.  For binary variables (maybe).  Also for matrix
equations.

### Graphics

For animation.

### Dynamics

For setting up and solving (numerically) dynamical system differential
equations.

### Averaging

To perform averaging of dynamical systems.


### Differential Geometry

For setting up and manipulating systems defined on manifolds.  Also has
basic stuff, like vectors and co-vectors, and defines the valid
operations between these constructs.

### Geometric Mechanics

Geometric mechanics here is mainly about working out the geometry of
dynamical systems on Lie Groups.  There is related code in the Mechanics
package, which looks at Lagrangian or Hamiltonian systems.


### Mechanics

In the sense of Lagrangian and Hamiltonian mechanics.  Mostly as would
be used to simulate and describe physical systems.


### Robotics

Stuff built up when taking Robotics course at Caltech (ME115).  Probably
not the best since I was learning Mathematica at the same time.

