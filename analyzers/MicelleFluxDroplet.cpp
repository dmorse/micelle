#ifndef MCMD_MICELLE_FLUX_DROPLET_CPP
#define MCMD_MICELLE_FLUX_DROPLET_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MicelleFluxDroplet.h"
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
   MicelleFluxDroplet::MicelleFluxDroplet(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      outputFile_(),
      dropletSpeciesId_(),
      surfactantSpeciesId_(),
      dropletAtomTypeId_(),
      surfactantAtomTypeId_(),
      cutoff_(),
      histMin_(),
      histMax_(),
      beadNumber_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("MicelleFlux"); }

   /// Read parameters from file, and allocate arrays.
   void MicelleFluxDroplet::readParameters(std::istream& in) 
   {   
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "dropletSpeciesId", dropletSpeciesId_);
      if (dropletSpeciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      read<int>(in, "surfactantSpeciesId", surfactantSpeciesId_);
      if (surfactantSpeciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (surfactantSpeciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      read<int>(in, "dropletAtomTypeId", dropletAtomTypeId_);
      if (dropletAtomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      read<int>(in, "surfactantAtomTypeId", surfactantAtomTypeId_);
      if (surfactantAtomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      if (surfactantAtomTypeId_ >= system().simulation().nAtomType()) {
         UTIL_THROW("nTypeId >= nAtomType");
      }

      read<double>(in, "cutoff", cutoff_);

      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      read<double>(in, "radius", radius_);

      read<int>(in, "beadNumber", beadNumber_);
      read<int>(in, "dropletLength", dropletLength_);
      isInitialized_ = true;
      /*
      // Initialize ClusterIdentifier
      std::cout << "YES";
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
      */
   }

   /*
   * Load state from an archive.
   */
   void MicelleFluxDroplet::loadParameters(Serializable::IArchive& ar)
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
   void MicelleFluxDroplet::save(Serializable::OArchive& ar)
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
   void MicelleFluxDroplet::setup() 
   {    
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      nSample_ = 0;
      identifier_.initialize(dropletSpeciesId_, dropletAtomTypeId_, cutoff_);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
      int surfactantSpeciesCapacity=system().simulation().species(surfactantSpeciesId_).capacity();
      int dropletSpeciesCapacity=system().simulation().species(dropletSpeciesId_).capacity();
      micelleFlux_.allocate(surfactantSpeciesCapacity);
      priorMicelleFlux_.allocate(surfactantSpeciesCapacity);
      inDropletDroplet_.allocate(dropletSpeciesCapacity);
      inDropletSurfactant_.allocate(surfactantSpeciesCapacity);
      for (int i = 0; i < surfactantSpeciesCapacity; i++) {
        priorMicelleFlux_[i]=-1;
      }
      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);
   }

   /*
   *  Calculate the micelle's center of mass
   */
   
   Vector MicelleFluxDroplet::comCalculator(DArray<int> micelleIds)
   { 
     Species* speciesPtr;
     speciesPtr = &(system().simulation().species(dropletSpeciesId_));
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
         for (int k = 0; k < dropletLength_; ++k) { 
           particleCount_ = particleCount_+1;
           r = system().molecule(dropletSpeciesId_, i).atom(k).position();
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
   void MicelleFluxDroplet::sample(long iStep) 
   {   
      if (isAtInterval(iStep)) {
      // Step one, determine which molecules make up the droplet among droplet forming molecules
         identifier_.identifyClusters();
         ++nSample_;
      
      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(dropletSpeciesId_)); 
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
         if (clusterSize > 10) {
         inDropletDroplet_[i]=1;
         }
         else {
         inDropletDroplet_[i]=0;  
         }
      }
      Vector r;
      // Determine the location of the droplet center of mass
      micelleCOM_=comCalculator(inDropletDroplet_);

      //Step two determine which surfactant molecules are in contact with the droplet? or `

      cellList_.setup(system().boundary(), cutoff_);
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      CellList::NeighborArray neighborArray;
      Atom* otherAtom;
      double cutoffSq = cutoff_*cutoff_;
      double rsq;
      bool isMoleculeOnDroplet;
      int i = 0;
      system().begin(dropletSpeciesId_,molIter);
      for( ; molIter.notEnd(); ++molIter) { 
        for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
          if (atomIter->typeId() == dropletAtomTypeId_ and inDropletDroplet_[i] == 1  ) {
            system().boundary().shift(atomIter->position());
            cellList_.addAtom(*atomIter);
          }
        }
        i += 1;
      }
      /// Cell list built, now cycle through the other species to see which are in contact
      system().begin(surfactantSpeciesId_,molIter);
      i = 0;
      for( ; molIter.notEnd(); ++molIter) {
        isMoleculeOnDroplet = false;
        for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
          if (atomIter->typeId() == surfactantAtomTypeId_) {
            cellList_.getNeighbors(atomIter->position(), neighborArray);
            for (int i = 0; i< neighborArray.size(); i++) {
              otherAtom = neighborArray[i];
              rsq = system().boundary().distanceSq(atomIter->position(),otherAtom->position());
              if (rsq < cutoffSq) {
                isMoleculeOnDroplet = true;
                break;
              }  
            }
            if (isMoleculeOnDroplet) {
              inDropletSurfactant_[i] = 1;
              break;
            } else {
              inDropletSurfactant_[i]=0;
            }
          }
        }
        i += 1;
      }

     
      // Step three determine which surfactant molecules are outside the outer radius
      speciesPtr = &(system().simulation().species(surfactantSpeciesId_)); 
      nMolecules = speciesPtr->capacity();
      double distance;
      for (int i = 0; i < nMolecules; i++) {
        r = system().molecule(surfactantSpeciesId_, i).atom(beadNumber_).position();
        distance = sqrt(system().boundary().distanceSq(r,micelleCOM_));
        //Current trying revers
        if (inDropletSurfactant_[i]==1) {
          micelleFlux_[i] = 0;
        }
        else if (distance>radius_) {
          micelleFlux_[i] = 1;
        }
        else {
          micelleFlux_[i] = priorMicelleFlux_[i];
        }
      }
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
   void MicelleFluxDroplet::output() 
   {  /*  
      outputFile_.close();
      // Write parameter file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Write histogram output
      // hist_.output(outputFile_);
      */
   }

}
#endif 
