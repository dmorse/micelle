#include "MicelleMcModule.h"

namespace McMd 
{

   using namespace Util;
   
   MicelleMcModule::MicelleMcModule(McSimulation& sim)
    : analyzerFactory_(sim, sim.system()),
      mcMoveFactory_(sim, sim.system())
   {
      sim.analyzerFactory().addSubfactory(analyzerFactory_);
      sim.mcMoveFactory().addSubfactory(mcMoveFactory_);
   }

}
