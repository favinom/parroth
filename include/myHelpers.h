#pragma once

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

typedef libMesh::TensorValue<double> RealTensorValue;
typedef libMesh::VectorValue<double> RealVectorValue;

void fillRowOfTensor(RealVectorValue const & r,int const comp,RealTensorValue & V);

