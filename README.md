# Skeletonization by propagation

This software is an implementation of the skeletonization by propagation method.

The propagated skeleton: a robust detail-preserving approach, Durix B., Chambon S., Leonard K., Mari J.-L. and Morin G., Discrete Geometry for Computer Imagery, 2019

**Note:** An improved skeletonization method is available [here](https://github.com/Ibujah/compactskel).

## Website:

http://durix.perso.enseeiht.fr/

## Needed libraries:

 * [Boost](http://www.boost.org/)
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (header only)
 * [OpenCV 3](http://opencv.org/)
 * [Nlopt](http://ab-initio.mit.edu/nlopt)

## Instructions

Install the needed libraries

```
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

The program is now in bin/

To use it on an example:

```
./soft_2dskeletonization --img ../ressources/rat.png --output
```

To compare the propagated skeleton with some pruning methods:

```
./soft_2dskeletonization --img ../ressources/rat.png --compare --output
```
