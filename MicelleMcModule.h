#ifndef MCMD_MICELLE_MC_MODULE_H
#define MCMD_MICELLE_MC_MODULE_H

// Custom factory classes
#include <mcMd/mcSimulation/McSimulation.h>
#include "analyzers/MicelleAnalyzerFactory.h"
#include "mcMoves/MicelleMcMoveFactory.h"

namespace McMd 
{

   using namespace Util;

   /**
   * Module for slip link simulation. 
   */
   class MicelleMcModule 
   {

   public:

      /**
      * Constructor.  
      *
      * \param sim parent McSimulation.
      */
      MicelleMcModule(McSimulation& sim);

   private:
   
      MicelleAnalyzerFactory analyzerFactory_;
      MicelleMcMoveFactory mcMoveFactory_;
   
   };

}
#endif
