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
1. Create and cd into a \textbf{build}
directory in the \textbf{shape-lme} package part of the tutorial.
2. Run \textbf{ccmake ../} and \textbf{make}
