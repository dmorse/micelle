#ifndef MCMD_RADIAL_COMPOSITION_CPP
#define MCMD_RADIAL_COMPOSTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RadialComposition.h"
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
#include <simp/boundary/Boundary.h>
#include <sstream>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   RadialComposition::RadialComposition(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("RadialComposition"); }

   /// Read parameters from file, and allocate arrays.
   void RadialComposition::readParameters(std::istream& in) 
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

      read<double>(in, "radiusMax", radius_);

      read<int>(in, "nBins", nBins_);
      // Initialize ClusterIdentifier
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      isInitialized_ = true;
      int speciesCapacity=system().simulation().species(speciesId_).capacity();
      InMicelle_.allocate(speciesCapacity);
   }

   /*
   * Load state from an archive.
   */
   void RadialComposition::loadParameters(Serializable::IArchive& ar)
   {
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
      
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      ar >> nSample_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void RadialComposition::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & nSample_;
   }

   /*
   * Clear accumulators.
   */
   void RadialComposition::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      nSample_ = 0;
      double min = 0.0;
      double max = radius_;
      nAtomType_ = system().simulation().nAtomType();
      accumulators_.allocate(nAtomType_);
      for (int i=0; i<nAtomType_;++i) {
         accumulators_[i].setParam(min,max,nBins_);
      }
      for (int i=0; i<nAtomType_;++i) {
        accumulators_[i].clear();
      }
   }

   /*
   *  Calculate the micelle's center of mass
   */
   
   Vector RadialComposition::comCalculator(DArray<int> micelleIds)
   { Species* speciesPtr;
     speciesPtr = &(system().simulation().species(speciesId_));
     int nMolecules = speciesPtr -> capacity();
     int beadsInChain = system().simulation().species(speciesId_).nAtom();
     int clusterSize = 0;
//     Vector centralMolecule;
     Vector r;
     Vector comTrack;
     for (int j=0; j<Dimension; ++j){
       comTrack[j] =0;
     }
     Vector centralAtom;
     Vector lengths = system().boundary().lengths();
     particleCount_=0;  
     for (int i = 0; i < nMolecules; ++i) {
       if (micelleIds[i] == 1)
         {  for (int k = 0; k < beadsInChain; ++k) {
            if (atomTypeId_==system().simulation().species(speciesId_).atomTypeId(k)) { 
            particleCount_ = particleCount_+1;
            r = system().molecule(speciesId_, i).atom(k).position();
            if (clusterSize == 0)
               {centralAtom = r;
               }
            clusterSize = clusterSize + 1;
            for (int j=0; j<Dimension; ++j) {
                if (std::abs(centralAtom[j]-r[j]) > lengths[j]/2) {
                    if ((r[j]-centralAtom[j]) > 0)
                      {comTrack[j] = comTrack[j]+r[j]-lengths[j];}
                      else 
                      {comTrack[j] = comTrack[j]+r[j]+lengths[j];}
                    } else {
                    comTrack[j]=comTrack[j]+r[j];
                }
              
            }
            }
           }
         }
      }  
         for (int j=0; j<Dimension; ++j) {
            comTrack[j]=comTrack[j]/particleCount_;
        }
        return comTrack;
   
   }


   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void RadialComposition::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {
         identifier_.identifyClusters();
         ++nSample_;
      Species* speciesPtr;
      System::ConstMoleculeIterator  molIter;
      Molecule::ConstAtomIterator  atomIter;
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
         if (clusterSize > 10) {
         InMicelle_[i]=1;
         }
         else {
         InMicelle_[i]=0;
         }
      }
      Vector r;
      // Calculate the micelle center of mass
      micelleCOM_=comCalculator(InMicelle_);
      Vector lengths = system().boundary().lengths();
      double distance;
      // Now for each atom in the system of the desginated type measure the distance from the COM and put it in the box it belongs in
      int nSpecies = system().simulation().nSpecies();
      int typeId;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        system().begin(iSpecies,molIter);
      for (;molIter.notEnd();++molIter) {
          molIter->begin(atomIter);
        for (;atomIter.notEnd(); ++atomIter) {
            distance = 0;
            r = atomIter->position();
            typeId = atomIter->typeId();
            for (int j = 0; j < Dimension; j++)
            { if (std::abs(micelleCOM_[j]-r[j]) > lengths[j]/2) {
                      if ((r[j]-micelleCOM_[j]) > 0)
                        distance = distance + pow(micelleCOM_[j]-r[j]+lengths[j],2);
                        else 
                        distance = distance + pow(micelleCOM_[j]-r[j]-lengths[j],2);
                      } else {
                      distance = distance + pow(micelleCOM_[j]-r[j],2);
                  }
            }
            distance = sqrt(distance);
            accumulators_[typeId].sample(distance);
      }
     }
    }
    } 
   }

   /*
   * Output results to file after simulation is completed.
   */
   void RadialComposition::output() 
   {  outputFile_.close();
      // Write parameter file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      // Write histogram output
      for (int i =0; i < nAtomType_; ++i) {
      accumulators_[i].output(outputFile_);
      outputFile_ << std::endl;
      }
      outputFile_.close();
   }

}
#endif 
