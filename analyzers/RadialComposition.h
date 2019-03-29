#ifndef MCMD_RADIAL_COMPOSITION_H
#define MCMD_RADIAL_COMPOSITION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>       // base class template
#include <mcMd/simulation/System.h>                   // class templ param
#include <mcMd/analyzers/system/ClusterIdentifier.h>  // member
#include <simp/boundary/Boundary.h>                   // member (typedef)
#include <util/accumulators/Distribution.h>           // member
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/containers/DArray.h>

namespace McMd
{

   using namespace Util;

   /**
   * Identify micelle clusters in polymeric systems.
   */
   class RadialComposition : public SystemAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      RadialComposition(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval         : sampling interval
      *   - string outputFileName   : base name for output file(s)
      *   - int    speciesId        : integer id for Species that makes up the micelle
      *   - int    coreId           : integer id for the type that makes up teh core 
      *   - double cutoff           : distance cutoff
      *   - int    nBins            : number of bins 
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /**
      * Identify the micelle, its  COM, and then radial composition around it
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
      
      /**
      *   Calculate the COM for the micelle
      *   
      *
      **/

      Vector comCalculator(DArray<int> micelleIds); 

      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// ClusterIdentifier
      ClusterIdentifier identifier_;
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      int  atomTypeId_;

      /// Distance cutoff
      double cutoff_;

      //Number of bins in the radial composition function
      int nBins_;
      /// distance outside the cluster required for a molecule to be considered outside the cluster 
      double radius_;

      /// Number of configurations dumped thus far (first dump is zero).
      long  nSample_;

      /// Has readParam been called?
      bool  isInitialized_;

      /// Is the species of interest mutable? 
      bool  isMutable_;
      /// If the species is mutable what subtype is of interest?
      int   speciesSubtype_;
      
      DArray<int> InMicelle_;
      //// Various stuff for position calcs
      // The micelle COM
      Vector micelleCOM_;
  
      // Number of atom types in the system
      int nAtomType_;
    
      DArray<Distribution> accumulators_; 
      int particleCount_;
      DArray<Vector> SurfactantPositions_;
      DArray<Vector> UnwrappedPositions_;
   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void RadialComposition::serialize(Archive& ar, const unsigned int version)
   {  
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & nSample_;
   }

}
#endif
