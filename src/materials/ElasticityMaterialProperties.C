//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticityMaterialProperties.h"
#include "myHelpers.h"

registerMooseObject("MooseApp", ElasticityMaterialProperties);

InputParameters
ElasticityMaterialProperties::validParams()
{

  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a property using a piecewise linear interpolation to define "
                             "its dependence on a variable");
  params.addRequiredParam<Real>("mu", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("lambda", "Scale factor to be applied to the ordinate values");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");

  return params;
}

ElasticityMaterialProperties::ElasticityMaterialProperties(const InputParameters & parameters) :
Material(parameters),
_dim(_mesh.dimension()),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),
_mu_prop(declareProperty<Real>("mu_prop")),
_lambda_prop(declareProperty<Real>("lambda_prop")),
_sigma_prop(declareProperty<RealTensorValue>("sigma_prop"))
{}

void
ElasticityMaterialProperties::computeQpProperties()
{
  _mu_prop[_qp]=_mu;
  _lambda_prop[_qp]=_lambda;

  _U.zero();
  fillRowOfTensor(_grad_disp_x[_qp],0,_U);
  fillRowOfTensor(_grad_disp_y[_qp],1,_U);
  fillRowOfTensor(_grad_disp_z[_qp],2,_U);

  computeStress(_U,_mu_prop[_qp],_lambda_prop[_qp],_sigma_prop[_qp]);
}
