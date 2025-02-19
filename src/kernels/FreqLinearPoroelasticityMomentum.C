//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FreqLinearPoroelasticityMomentum.h"
#include "myHelpers.h"

registerMooseObject("parrothApp", FreqLinearPoroelasticityMomentum);

InputParameters
FreqLinearPoroelasticityMomentum::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<bool>("real", "Component");
  params.addRequiredParam<int>("component", "Component");

  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("p", "pressure");

  return params;
}

FreqLinearPoroelasticityMomentum::FreqLinearPoroelasticityMomentum(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_is_real(getParam<bool>("real")),
_component(getParam<int>("component")),
_mu(getMaterialProperty<Real>("mu_prop")),
_lambda(getMaterialProperty<Real>("lambda_prop")),
_alpha(getMaterialProperty<Real>("alpha_prop")),
_sigma_r(getMaterialProperty<RealTensorValue>("sigma_r_prop")),
_sigma_i(getMaterialProperty<RealTensorValue>("sigma_i_prop")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p"))
{
	_I.zero();
	for(int i=0;i<_dim;++i)
		_I(i,i)=1.0;
}

Real
FreqLinearPoroelasticityMomentum::computeQpResidual()
{
	_V.zero();
	fillRowOfTensor(_grad_test[_i][_qp],_component,_V);
	if (_is_real)
		_sigmaH=_sigma_r[_qp];
	else
		_sigmaH=_sigma_i[_qp];

	return _sigmaH.contract(_V);
}


Real
FreqLinearPoroelasticityMomentum::computeQpJacobian()
{
	_H.zero();
	_V.zero();
	fillRowOfTensor(_grad_phi[_j][_qp],_component,_H);
	fillRowOfTensor(_grad_test[_i][_qp],_component,_V);

	computeStress(_H,_mu[_qp],_lambda[_qp],_sigmaH);

  return _sigmaH.contract(_V);
}


Real
FreqLinearPoroelasticityMomentum::computeQpOffDiagJacobian(unsigned int jvar) 
{
int componentH=-1;

if(jvar == _id_x)
	componentH = 0;

if(jvar == _id_y)
	componentH = 1;

if(jvar == _id_z)
	componentH = 2;

if(jvar == _id_p)
	componentH = 3;

_V.zero();
fillRowOfTensor(_grad_test[_i][_qp],_component,_V);

_sigmaH.zero();

if (0<=componentH & componentH<=2)
{
	_H.zero();
	fillRowOfTensor(_grad_phi[_j][_qp],componentH,_H);
	computeStress(_H,_mu[_qp],_lambda[_qp],_sigmaH);
}
if (componentH==3)
{
	_sigmaH=-_alpha[_qp]*_phi[_j][_qp]*_I;
}

return _sigmaH.contract(_V);

}


