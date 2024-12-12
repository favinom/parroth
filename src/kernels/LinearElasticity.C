//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearElasticity.h"

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
_component(getParam<int>("component")),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999)
 {
 std::cout << _dim << std::endl;
 }

Real
LinearElasticity::computeQpResidual()
{
RealTensorValue U;
RealTensorValue eps;
RealTensorValue I;
RealTensorValue sigma;
RealTensorValue V;
Real treps;
U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);
eps = (U+U.transpose())/2.0;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	I(i,j) = (i==j);
	V(i,j) = 0;
	}
}
if(_dim == 2){
	I(2,2) = 0.0;
}

treps = eps.tr();

sigma = 2*_mu[_qp]*eps + _lambda[_qp]*treps*I;

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

  return sigma.contract(V);
}


Real
LinearElasticity::computeQpJacobian()
{
RealTensorValue H;
RealTensorValue epsH;
RealTensorValue I;
RealTensorValue V;
RealTensorValue sigmaH;
Real trepsH;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	I(i,j) = (i==j);
	H(i,j) = 0;
	V(i,j) = 0;
	}
}

H(_component,0) = _grad_phi[_j][_qp](0); H(_component,1) = _grad_phi[_j][_qp](1); H(_component,2) = _grad_phi[_j][_qp](2);

if(_dim == 2){
	I(2,2) = 0.0;
	H(_component,2) = 0.0;
}

epsH = (H+H.transpose())/2.0;

trepsH = epsH.tr();

sigmaH = 2*_mu[_qp]*epsH + _lambda[_qp]*trepsH*I;

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

  return sigmaH.contract(V);
}


Real
LinearElasticity::computeQpOffDiagJacobian(unsigned int jvar) 
{
RealTensorValue H;
RealTensorValue epsH;
RealTensorValue I;
RealTensorValue V;
RealTensorValue sigmaH;
Real trepsH;
int componentH;

if(jvar == _id_x){
	componentH = 0;
}

if(jvar == _id_y){
	componentH = 1;
}

if(jvar == _id_z){
	componentH = 2;
}

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	I(i,j) = (i==j);
	H(i,j) = 0;
	V(i,j) = 0;
	}
}
if(_dim == 2){
	I(2,2) = 0.0;
}

H(componentH,0) = _grad_phi[_j][_qp](0); H(componentH,1) = _grad_phi[_j][_qp](1); H(componentH,2) = _grad_phi[_j][_qp](2);


epsH = (H+H.transpose())/2.0;

trepsH = epsH.tr();

sigmaH = 2*_mu[_qp]*epsH + _lambda[_qp]*trepsH*I;

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

return sigmaH.contract(V);

}


