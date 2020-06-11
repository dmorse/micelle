#ifndef MCMD_SINGLE_MICELLE_UMBRELLA_SAMPLING_MOVE_H
#define MCMD_SINGLE_MICELLE_UMBRELLA_SAMPLING_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class
#include <util/containers/DArray.h>
#include <util/containers/Pair.h>
#include <mcMd/chemistry/Molecule.h> //testing
#include <mcMd/modules/micelle/mcMoves/micelle/ClusterIdentifierMC.h>
#include <util/accumulators/IntDistribution.h> 

namespace McMd
{

   using namespace Util;
   class SpeciesMutator;
   class LinearSG;
   class McSystem;

   /**
   * A move that changes the type of a HomopolymerSG molecule.
   *
   * \ingroup McMd_McMove_Module
   */
   class SingleMicelleUmbrellaSamplingMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      SingleMicelleUmbrellaSamplingMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
    
      Molecule& randomSGMolecule(int speciesId, int typeId, int flipType, bool atBound);
      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
      
      virtual void output();

  protected:
      /// Integer index for molecular species.
      int speciesId_;

      int uLimit_;
      int flipsAttempted_, flipsAccepted_, MicelleRejections_;

      IntDistribution hist_;

      int lLimit_;
 
      int flipSubType_;

      int flipsPerMove_;
   
      int aggregateCutoff_;

      DArray<int> oldIdentities_;
 
      DArray<double> weights_;

      DArray<int> stateCount_;

      /// Pointer to instance of HomopolymerSG.
      LinearSG* speciesPtr_;
      SpeciesMutator* mutatorPtr_;

      std::string outputFileName_;
      std::string initialWeightFileName_;
      std::string clusterHistName_;
      DArray<int> steps_;
      DArray<double> weightTrack_;
      int stepCount_;
      int flipType_;
      std::ofstream outputFile_;
      std::ofstream clusterFile_;
      //Total number of semigrand molecules possible
      int capacity_;
      ClusterIdentifierMC identifier_;
      double cutoff_;
      int atomTypeId_;
      int histMin_;
      int histMax_;
   };

}      
#endif
