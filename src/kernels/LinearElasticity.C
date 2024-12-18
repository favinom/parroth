//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearElasticity.h"
#include "myHelpers.h"

registerMooseObject("parrothApp", LinearElasticity);

InputParameters
LinearElasticity::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<int>("component", "Component");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  //params.addCoupledVar("density",
  //                     "The name of the temperature variable used in the "
  //                     "ComputeThermalExpansionEigenstrain.  (Not required for "
  //                     "simulations without temperature coupling.)");
  return params;
}

LinearElasticity::LinearElasticity(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_mu(getMaterialProperty<Real>("mu_prop")),
_lambda(getMaterialProperty<Real>("lambda_prop")),
_sigma(getMaterialProperty<RealTensorValue>("sigma_prop")),
_component(getParam<int>("component")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999)
{
	_I.zero();
	for(int i=0;i<3;++i)
		_I(i,i)=1.0;
}

Real
LinearElasticity::computeQpResidual()
{
	_V.zero();
	fillRowOfTensor(_grad_test[_i][_qp],_component,_V);
	return _sigma[_qp].contract(_V);
}


Real
LinearElasticity::computeQpJacobian()
{
	_H.zero();
	_V.zero();
	fillRowOfTensor(_grad_phi[_j][_qp],_component,_H);
	fillRowOfTensor(_grad_test[_i][_qp],_component,_V);

	computeStress(_H,_mu[_qp],_lambda[_qp],_sigmaH);

  return _sigmaH.contract(_V);
}


Real
LinearElasticity::computeQpOffDiagJacobian(unsigned int jvar) 
{
int componentH;

if(jvar == _id_x)
	componentH = 0;

if(jvar == _id_y)
	componentH = 1;

if(jvar == _id_z)
	componentH = 2;

_H.zero();
_V.zero();
fillRowOfTensor(_grad_phi[_j][_qp],componentH,_H);
fillRowOfTensor(_grad_test[_i][_qp],_component,_V);

computeStress(_H,_mu[_qp],_lambda[_qp],_sigmaH);
return _sigmaH.contract(_V);

}


