# PHD_Shape
A private repository for coordinating efforts with the Utah/NYU groups

### Deformetrica Build Instructions
1. Build ITK
  * Turn on 'ITK_WRAP' 
  * Turn on 'FFTWD'
  * Turn on 'FFTWF'
2. Build VTK
3. Build Armadillo
  * http://arma.sourceforge.net/
4. Build Deformetrica against above three.

### ShapeLME Build Instructions 
Dependencies: Armadillo (already built for Deformetrica)

* mkdir shape-lme/bin and cd into it.
* Run ccmake ../ and follow instructions.
* Run make
