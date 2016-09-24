/*==============================================================================
  File: primitives.hxx

  Typedefs for primitives 

==============================================================================*/

#ifndef __primitives_hxx
#define __primitives_hxx

#include <vector>
#include <armadillo>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <memory>

using namespace arma;
using namespace std;

// basic primitives
typedef double ScalarType;

typedef Col<ScalarType> ScalarListType;
typedef Col<unsigned int> IntListType;
typedef Mat<ScalarType> MatrixType;
typedef Col<ScalarType> VectorType;
typedef Row<ScalarType> RowVectorType;

typedef vector<unsigned int> TimePointsType;

#endif
