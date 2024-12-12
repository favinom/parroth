//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "parrothTestApp.h"
#include "parrothApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
parrothTestApp::validParams()
{
  InputParameters params = parrothApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

parrothTestApp::parrothTestApp(InputParameters parameters) : MooseApp(parameters)
{
  parrothTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

parrothTestApp::~parrothTestApp() {}

void
parrothTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  parrothApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"parrothTestApp"});
    Registry::registerActionsTo(af, {"parrothTestApp"});
  }
}

void
parrothTestApp::registerApps()
{
  registerApp(parrothApp);
  registerApp(parrothTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
parrothTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrothTestApp::registerAll(f, af, s);
}
extern "C" void
parrothTestApp__registerApps()
{
  parrothTestApp::registerApps();
}
