//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FreqPoroelasticityMaterialProperties.h"
#include "myHelpers.h"

registerMooseObject("MooseApp", FreqPoroelasticityMaterialProperties);

InputParameters
FreqPoroelasticityMaterialProperties::validParams()
{

  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a property using a piecewise linear interpolation to define "
                             "its dependence on a variable");
  params.addRequiredParam<Real>("mu", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("lambda", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("alpha", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("m", "Scale factor to be applied to the ordinate values");
  params.addRequiredParam<Real>("k", "Scale factor to be applied to the ordinate values");

  params.addRequiredCoupledVar("disp_r_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_i_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_r_y", "Displacement along y");
  params.addRequiredCoupledVar("disp_i_y", "Displacement along y");
  params.addCoupledVar("disp_r_z","Displacement along z");
  params.addCoupledVar("disp_i_z","Displacement along z");
  params.addRequiredCoupledVar("p_r", "Displacement along y");
  params.addRequiredCoupledVar("p_i", "Displacement along y");


  return params;
}

FreqPoroelasticityMaterialProperties::FreqPoroelasticityMaterialProperties(const InputParameters & parameters) :
Material(parameters),
_dim(_mesh.dimension()),
_grad_disp_r_x(coupledGradient("disp_r_x")),
_grad_disp_i_x(coupledGradient("disp_i_x")),
_grad_disp_r_y(coupledGradient("disp_r_y")),
_grad_disp_i_y(coupledGradient("disp_i_y")),
_grad_disp_r_z(_dim>2 ? coupledGradient("disp_r_z") : _grad_zero),
_grad_disp_i_z(_dim>2 ? coupledGradient("disp_i_z") : _grad_zero),
_p_r(coupledValue("p_r") ),
_p_i(coupledValue("p_i") ),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda")),
_alpha(getParam<Real>("alpha")),
_m(getParam<Real>("m")),
_k(getParam<Real>("k")),
_mu_prop(declareProperty<Real>("mu_prop")),
_lambda_prop(declareProperty<Real>("lambda_prop")),
_alpha_prop(declareProperty<Real>("alpha_prop")),
_m_prop(declareProperty<Real>("m_prop")),
_k_prop(declareProperty<Real>("k_prop")),
_eps_r_prop(declareProperty<RealTensorValue>("eps_r_prop")),
_eps_i_prop(declareProperty<RealTensorValue>("eps_i_prop")),
_sigma_r_prop(declareProperty<RealTensorValue>("sigma_r_prop")),
_sigma_i_prop(declareProperty<RealTensorValue>("sigma_i_prop"))
{
  _I.zero();
  for (int i=0; i<_dim; ++i)
    _I(i,i)=1.0;
}

void
FreqPoroelasticityMaterialProperties::computeQpProperties()
{
  _mu_prop[_qp]=_mu;
  _lambda_prop[_qp]=_lambda;
  _alpha_prop[_qp]=_alpha;
  _m_prop[_qp]=_m;
  _k_prop[_qp]=_k;

  _U_r.zero();
  _U_i.zero();
  fillRowOfTensor(_grad_disp_r_x[_qp],0,_U_r);
  fillRowOfTensor(_grad_disp_r_y[_qp],1,_U_r);
  fillRowOfTensor(_grad_disp_r_z[_qp],2,_U_r);

  fillRowOfTensor(_grad_disp_i_x[_qp],0,_U_i);
  fillRowOfTensor(_grad_disp_i_y[_qp],1,_U_i);
  fillRowOfTensor(_grad_disp_i_z[_qp],2,_U_i);

  _eps_r_prop[_qp]=0.5*(_U_r+_U_r.transpose());
  _eps_i_prop[_qp]=0.5*(_U_i+_U_i.transpose());

  computeStress(_U_r,_mu_prop[_qp],_lambda_prop[_qp],_sigma_r_prop[_qp]);
  computeStress(_U_i,_mu_prop[_qp],_lambda_prop[_qp],_sigma_i_prop[_qp]);

  _sigma_r_prop[_qp]=_sigma_r_prop[_qp]-_alpha_prop[_qp]*_p_r[_qp]*_I;
  _sigma_i_prop[_qp]=_sigma_i_prop[_qp]-_alpha_prop[_qp]*_p_i[_qp]*_I;

}
