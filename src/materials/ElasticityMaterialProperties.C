//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticityMaterialProperties.h"

registerMooseObject("MooseApp", ElasticityMaterialProperties);

InputParameters
ElasticityMaterialProperties::validParams()
{

  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a property using a piecewise linear interpolation to define "
                             "its dependence on a variable");
  params.addRequiredParam<Real>("mu", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("lambda", "Scale factor to be applied to the ordinate values");
  return params;
}

ElasticityMaterialProperties::ElasticityMaterialProperties(const InputParameters & parameters) :
Material(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),
_mu_prop(declareProperty<Real>("mu_prop")),
_lambda_prop(declareProperty<Real>("lambda_prop"))
{}

void
ElasticityMaterialProperties::computeQpProperties()
{
  _mu_prop[_qp]=_mu;
  _lambda_prop[_qp]=_lambda;
}
