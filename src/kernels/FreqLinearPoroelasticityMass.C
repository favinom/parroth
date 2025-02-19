//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FreqLinearPoroelasticityMass.h"
#include "myHelpers.h"

registerMooseObject("parrothApp", FreqLinearPoroelasticityMass);

InputParameters
FreqLinearPoroelasticityMass::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<bool>("real", "Component");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("other_p", "pressure");

  return params;
}

FreqLinearPoroelasticityMass::FreqLinearPoroelasticityMass(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_is_real(getParam<bool>("real")),
_alpha(getMaterialProperty<Real>("alpha_prop")),
_m(getMaterialProperty<Real>("m_prop")),
_k(getMaterialProperty<Real>("k_prop")),
_eps_r(getMaterialProperty<RealTensorValue>("eps_r_prop")),
_eps_i(getMaterialProperty<RealTensorValue>("eps_i_prop")),
_other_p(coupledValue("other_p")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_other_p(coupled("other_p"))
{}

Real
FreqLinearPoroelasticityMass::computeQpResidual()
{
	Real omega=2.0*pi*std::pow(10.0,_t);
	Real ret=-(_k[_qp]/omega)*(_grad_u[_qp]*_grad_test[_i][_qp]);
	//Real ret2;
	if (_is_real)
		ret=ret+(_m[_qp]*_other_p[_qp]+_alpha[_qp]*_eps_i[_qp].tr())*_test[_i][_qp];
	else
		ret=ret-(_m[_qp]*_other_p[_qp]+_alpha[_qp]*_eps_r[_qp].tr())*_test[_i][_qp];


	return ret;
}


Real
FreqLinearPoroelasticityMass::computeQpJacobian()
{
	Real omega=2.0*pi*std::pow(10.0,_t);
	return -(_k[_qp]/omega)*(_grad_phi[_j][_qp]*_grad_test[_i][_qp]);
}


Real
FreqLinearPoroelasticityMass::computeQpOffDiagJacobian(unsigned int jvar) 
{
	Real ret=0.0;

	int componentH=-1;

	if(jvar == _id_x)
		componentH = 0;

	if(jvar == _id_y)
		componentH = 1;

	if(jvar == _id_z)
		componentH = 2;

	if(jvar == _id_other_p)
		componentH = 3;
	
	if (0<=componentH & componentH<=2)
		ret=_alpha[_qp]*_grad_phi[_j][_qp](componentH)*_test[_i][_qp];

	if (componentH==3)
		ret=_m[_qp]*_phi[_j][_qp]*_test[_i][_qp];

	if (_is_real==0)
		ret=-ret;

	return ret;

}

