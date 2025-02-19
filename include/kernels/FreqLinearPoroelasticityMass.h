//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class FreqLinearPoroelasticityMass : public Kernel
{
public:
  static InputParameters validParams();

  FreqLinearPoroelasticityMass(const InputParameters & parameters);
protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  int  const _dim;
  bool const _is_real;

  RealTensorValue _H;
  //RealTensorValue _V;
  //RealTensorValue _sigmaH;

  MaterialProperty<Real> const & _alpha;
  MaterialProperty<Real> const & _m;
  MaterialProperty<Real> const & _k;
  
  MaterialProperty<RealTensorValue> const & _eps_r;
  MaterialProperty<RealTensorValue> const & _eps_i;

  VariableValue const & _other_p;

  int const _id_x;
  int const _id_y;
  int const _id_z;
  int const _id_other_p;
  
};
