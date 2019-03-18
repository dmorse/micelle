#ifndef MICELLE_MC_MOVE_FACTORY_CPP
#define MICELLE_MC_MOVE_FACTORY_CPP

#include "MicelleMcMoveFactory.h" 
#include <mcMd/mcSimulation/McSystem.h> 

// Include headers for user defined McMoves
#include "micelle/SingleMicelleHybridMove.h"
#include "micelle/SingleMicelleUmbrellaSamplingMove.h"
#include "micelle/RgUmbrellaSamplingMove.h"

#include "semigrand/WangLandauSemigrandMove.h"
#include "semigrand/UmbrellaSamplingSemigrandMove.h"

namespace McMd
{


   /*
   * Constructor.
   */
   MicelleMcMoveFactory::MicelleMcMoveFactory(McSimulation& simulation, 
                                                McSystem& system)
    : McMoveFactory(simulation, system)
   {}

   /* 
   * Return a pointer to a new instance of className.
   */
   McMove* MicelleMcMoveFactory::factory(const std::string &className) const
   {
      McMove* ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "SingleMicelleHybridMove") {
         ptr = new SingleMicelleHybridMove(system());
      } else     
      if (className == "SingleMicelleUmbrellaSamplingMove") {
         ptr = new SingleMicelleUmbrellaSamplingMove(system());
      } else
      if (className == "RgUmbrellaSamplingMove") {
        ptr = new RgUmbrellaSamplingMove(system());
      } else          
      if (className == "WangLandauSemigrandMove") {
        ptr = new WangLandauSemigrandMove(system());
      } else          
      if (className == "UmbrellaSamplingSemigrandMove") {
        ptr = new UmbrellaSamplingSemigrandMove(system());
      }          
      
      #if 0     
      // If not a user-defined class, try the standard factory 
      if (!ptr) {
         ptr = McMoveFactory::factory(className);
      }
      #endif

      return ptr;
   }

}

#endif
