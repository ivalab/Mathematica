IvaMatica
===========

A collection of IVALab's Mathematica homebrew packages, mostly related to
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

**Others.**
There are other folk who have implemented stuff like this in
Mathematica.  For example, the 
[http://www.cds.caltech.edu/~murray/mlswiki/index.php/Software](software
associated) to the Murray, Li, Sastry book *A Mathematical Introduction
to Robotic Manipulation* which is a pretty standard geometric approach
to the mechanics of manipulation.  Another place might be the 
[http://motion.me.ucsb.edu/book-gcms/Mma/](packages) from the Bullo and
Lewis book *Geometric Control of Mechanical Systems.*  Both of these
examples have different perspectives in relation to the codebase here.
The former is more concerned with manipulators and manipulation, while
the latter is more concerned with Bullo and Lewis' conception of
mechanical systems, which tends to focus on the Riemmanian structure.
There is some overlap with regards to connection forms.
It does have some pretty good documentation and use cases.  If I can get
my library cleaned up and modernized, then maybe I can start to
integrate the complementary aspects from these other two libraries.

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
equations.  Includes also a package for performing averaging of
dynamical systems.

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

