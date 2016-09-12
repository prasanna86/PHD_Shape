# PHD_Shape
A private repository for coordinating efforts with the Utah/NYU groups

### Deformetrica Build Instruction
1. Build ITK
  * Turn on 'ITK_WRAP' 
  * Turn on 'FFTWD'
  * Turn on 'FFTWF'
2. Build VTK
3. Build Armadillo
  * http://arma.sourceforge.net/
4. Build Deformetrica against above three.

### ShapeLME Build Instruction 
Dependencies: Armadillo (already built for Deformetrica)
1. mkdir shape-lme/bin and cd into it.
2. Run ccmake ../ and follow instructions.
3. Run make
