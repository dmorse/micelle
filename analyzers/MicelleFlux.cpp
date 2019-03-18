#ifndef MCMD_MICELLE_FLUX_CPP
#define MCMD_MICELLE_FLUX_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MicelleFlux.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   MicelleFlux::MicelleFlux(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      histMin_(),
      histMax_(),
      beadNumber_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("MicelleFlux"); }

   /// Read parameters from file, and allocate arrays.
   void MicelleFlux::readParameters(std::istream& in) 
   {   
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(speciesId_)); 
      isMutable_ = (speciesPtr->isMutable());
      if (isMutable_) {
         read<int>(in, "speciesSubtype", speciesSubtype_);
      }
      read<int>(in, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      if (atomTypeId_ >= system().simulation().nAtomType()) {
         UTIL_THROW("nTypeId >= nAtomType");
      }

      read<double>(in, "cutoff", cutoff_);

      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      read<double>(in, "radius", radius_);

      read<int>(in, "beadNumber", beadNumber_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void MicelleFlux::loadParameters(Serializable::IArchive& ar)
   {  /*
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(speciesId_)); 
      isMutable_ = (speciesPtr->isMutable());
      if (isMutable_) {
         loadParameter<int>(ar, "speciesSubtype", speciesSubtype_);
      }

      loadParameter<int>(ar, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      loadParameter<double>(ar, "radius", radius_);
      
      loadParameter<int>(ar, "beadNumber", beadNumber_);
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      ar >> nSample_;

      isInitialized_ = true;
      */
   }

   /*
   * Save state to archive.
   */
   void MicelleFlux::save(Serializable::OArchive& ar)
   {  /*
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & nSample_;
      */
   }

   /*
   * Clear accumulators.
   */
   void MicelleFlux::setup() 
   {    
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      nSample_ = 0;
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
      int speciesCapacity=system().simulation().species(speciesId_).capacity();
      micelleFlux_.allocate(speciesCapacity);
      priorMicelleFlux_.allocate(speciesCapacity);
      InMicelle_.allocate(speciesCapacity);
      for (int i = 0; i < speciesCapacity; i++) {
        priorMicelleFlux_[i]=-1;
      }
   }

   /*
   *  Calculate the micelle's center of mass; this was written before individual clusters could do this
   */
   
   Vector MicelleFlux::comCalculator(DArray<int> micelleIds)
   { 
     Species* speciesPtr;
     speciesPtr = &(system().simulation().species(speciesId_));
     int nMolecules = speciesPtr -> capacity();
     int clusterSize = 0;
//     Vector centralMolecule;
     Vector r;
     Vector comTrack;
     comTrack.zero();
     Vector centralAtomPosition;
     Vector dr;
     particleCount_=0;  
     for (int i = 0; i < nMolecules; ++i) {
       if (micelleIds[i] == 1) {
         for (int k = 0; k < beadNumber_; ++k) { 
           particleCount_ = particleCount_+1;
           r = system().molecule(speciesId_, i).atom(k).position();
           if (clusterSize == 0) {
             centralAtomPosition = r;
             clusterSize = clusterSize + 1;
           }
           system().boundary().distanceSq(centralAtomPosition,r,dr);
           comTrack -= dr;
         }
       }
     } 
     for (int j=0; j<Dimension; ++j) {
       comTrack[j]=comTrack[j]/particleCount_;
     }
     comTrack+=centralAtomPosition;
     return comTrack;
     
   }


   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void MicelleFlux::sample(long iStep) 
   {   
      if (isAtInterval(iStep)) {
         identifier_.identifyClusters();
         ++nSample_;
      
      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(speciesId_)); 
      int nMolecules = speciesPtr->capacity();
      ClusterLink* LinkPtr;      
      // Write the time step to the data file
      // Cycle through all the surfactant molecules; if they are in a micelle set status = 1, otherwise set status = 0;
      for (int i = 0; i < nMolecules; i++) {
         bool inCluster = false;
         // Get the molecule link
         LinkPtr = &(identifier_.link(i));
         // Use the link to get the cluster id
         int clusterId = LinkPtr->clusterId();
         // Determine the cluster size
         int clusterSize = identifier_.cluster(clusterId).size();
         if (clusterSize > 15) {
         InMicelle_[i]=1;
         }
         else {
         InMicelle_[i]=0;
         
         }
      }
      Vector r;
      // Calculate micelle COM
      micelleCOM_=comCalculator(InMicelle_);
      // Now with the COM calculate each chain's distance from the COM; 
      // If not in the micelle, test if it is far enough away to be considered out
      // otherwise label it with teh same label as it had prior
      double distance;
      for (int i = 0; i < nMolecules; i++) {
        r = system().molecule(speciesId_, i).atom(beadNumber_).position();
        distance = sqrt(system().boundary().distanceSq(r,micelleCOM_));
        if (InMicelle_[i]==1) {
          micelleFlux_[i] = 0;
        }
        else if (distance>radius_) {
          micelleFlux_[i] = 1;
        }
        else {
          micelleFlux_[i] = priorMicelleFlux_[i];
        }
      }
      // Write all the states of chains
      for (int i = 0; i < nMolecules; i++) {
        outputFile_ << micelleFlux_[i] << "  ";
      if (i == nMolecules -1)
        outputFile_ << "\n";
      }
      priorMicelleFlux_=micelleFlux_;   
         
    }
    
     
   }

   /*
   * Output results to file after simulation is completed.
   */
   void MicelleFlux::output() 
   {    
      outputFile_.close();
   }

}
#endif 
