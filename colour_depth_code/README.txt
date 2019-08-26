Color-filtered aperture-based depth estimation code

        from

Yosuke Bando, Bing-Yu Chen and Tomoyuki Nishita.
"Extracting Depth and Matte using a Color-Filtered Aperture,"
ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia 2008),
Vol. 27, No. 5, Article 134 (9 pages), 2008.

        and

United States Patent Application 20090284627.


                                  Nov. 8, 2012
                                  Yosuke Bando
                                  ybando@nis-lab.is.s.u-tokyo.ac.jp
                                  yosuke1.bando@toshiba.co.jp


0. Preface

The code is provided without any warranty, and it may be
used for non-commercial purposes for free.



1. To build

The code depends on the MRF minimization library.
Download "MRFx.x.zip" (x.x is a version number)
from http://vision.middlebury.edu/MRF/code/,
and unzip it to a folder of your choice (say ./MRFx.x/).

Change the following two lines in mrf.h in the MRF library
    typedef int EnergyVal;        /* The total energy of a labeling */
    typedef int CostVal;          /* costs of individual terms of the energy */
to
    typedef float EnergyVal;        /* The total energy of a labeling */
    typedef float CostVal;          /* costs of individual terms of the energy */

and then build the library by typing "make" in the MRF folder.

Now go back to the folder containing this README,
and change the following line in Makefile appropriately
so that it points to the MRF folder.

    MRF       = ./MRFx.x

Typing "make" will create an executable named "main(.exe)".



2. For test run

% ./run.sh

Three images (.pgm/ppm) will be created.
- depth_local.ppm  : local depth estimation result (blue means no value)
- depth_global.pgm : depth map after global optimization via graph-cuts
- restored.ppm     : color misalignment-corrected image based on depth

Color correction here is just a supplementary output for
illustrating the validity of the depth estimation result,
and often contains errors around foreground object silhouette.
For better correction, we apply the operation separately to
foreground and background parts as described in the last
paragraph of Sec. 6 in the paper.



3. Description of the files

+ images/*.ppm : sample photos captured with color-filtered aperture camera
+ *.{h,c}      : source code of the depth estimation
+ Makefile     : Makefile for Linux/Cygwin
+ run.sh       : sample command for running the executable
+ README.txt   : this file



4. Command line arguments

main in.ppm size disp_min disp_max [options]
-d int  : downsample size
-s float: smoothness weight

Mandatory arguments
(1) in.ppm   : input image
(2) size     : local window size (integer in [1, 7], 7 was used in the paper)
               actual window size will be (size * 2 + 1) x (size * 2 + 1)
(3) disp_min : minimum disparity ( -5 in the paper)
(4) disp_max : maximum disparity (+10 in the paper)

Optional arguments
-d int
    This option specifies an image downsampling rate for graph-cut
    optimization, as a simple way to reduce computation time.
    For example, specifying 2 will produce a depth map with half
    the resolution both in width and height. Local depth estimation
    will still be performed in full resolution. Default is 1.

-s float
    This option controls the weight of the smoothness energy term
    in the graph-cut optimization. Default is 1.0. A larger value
    will lead to a smoother result with less noise but will also
    lead to loss of details.

