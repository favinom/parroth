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
class LinearElasticity : public Kernel
{
public:
  static InputParameters validParams();

  LinearElasticity(const InputParameters & parameters);
protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


	int const _dim;

  RealTensorValue _I;

	MaterialProperty<Real> const & _mu;
  MaterialProperty<Real> const & _lambda;


	int const _component;
    VariableGradient const & _grad_disp_x;
    VariableGradient const & _grad_disp_y;
    VariableGradient const & _grad_disp_z;
    int const _id_x;
    int const _id_y;
    int const _id_z;
  
};


