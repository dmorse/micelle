#include "MicelleAnalyzerFactory.h"  
#include <mcMd/mcSimulation/McSimulation.h>  
#include <mcMd/mcSimulation/McSystem.h>  

// Include headers for any user defined Analyzers for MC simulations
#include "ClusterHistogramInd.h"
#include "ClusterHistogramSimple.h"
#include "InterfacialLoading.h"
#include "MicelleFlux.h"
#include "MicelleFluxDroplet.h"
#include "RadialComposition.h"

namespace McMd
{

   /* 
   * Return a pointer to a new instance of className.
   */
   Analyzer* MicelleAnalyzerFactory::factory(const std::string &className) const
   {
      Analyzer* ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "ClusterHistogramInd") {
         ptr = new ClusterHistogramInd(system());
      }             
      else if (className == "ClusterHistogramSimple") {
         ptr = new ClusterHistogramSimple(system());
      }
      else if (className == "InterfacialLoading") {
         ptr = new InterfacialLoading(system());
      }      
      else if (className == "MicelleFlux") {
         ptr = new MicelleFlux(system());
      }  
      else if (className == "MicelleFluxDroplet") {
         ptr = new MicelleFluxDroplet(system());
      }         
      else if (className == "RadialComposition") {
         ptr = new RadialComposition(system());
      }

      #if 0
      // If not a user-defined class, try the standard factory 
      if (!ptr) {
         ptr = McAnalyzerFactory::factory(className);
      }
      #endif

      return ptr;
   }

}
