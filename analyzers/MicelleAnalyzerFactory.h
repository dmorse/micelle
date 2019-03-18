#ifndef MICELLE_MC_ANALYZER_FACTORY_H
#define MICELLE_MC_ANALYZER_FACTORY_H

#include <mcMd/analyzers/system/SystemAnalyzerFactory.h>

namespace McMd
{

   /**
   * Custom AnalyzerFactory for an McSimulation
   */
   class MicelleAnalyzerFactory : public SystemAnalyzerFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      MicelleAnalyzerFactory(McSimulation& simulation, McSystem& system)
       : SystemAnalyzerFactory(simulation, system)
      {}

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

   };

}
#endif
