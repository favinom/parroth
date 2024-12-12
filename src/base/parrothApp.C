#include "parrothApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
parrothApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

parrothApp::parrothApp(InputParameters parameters) : MooseApp(parameters)
{
  parrothApp::registerAll(_factory, _action_factory, _syntax);
}

parrothApp::~parrothApp() {}

void 
parrothApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<parrothApp>(f, af, s);
  Registry::registerObjectsTo(f, {"parrothApp"});
  Registry::registerActionsTo(af, {"parrothApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
parrothApp::registerApps()
{
  registerApp(parrothApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
parrothApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrothApp::registerAll(f, af, s);
}
extern "C" void
parrothApp__registerApps()
{
  parrothApp::registerApps();
}
