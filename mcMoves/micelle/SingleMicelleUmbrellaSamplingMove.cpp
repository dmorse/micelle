/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SingleMicelleUmbrellaSamplingMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/LinearSG.h>
#include <mcMd/chemistry/Molecule.h> //testing
#include <util/global.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   SingleMicelleUmbrellaSamplingMove::SingleMicelleUmbrellaSamplingMove(McSystem& system) : 
      SystemMove(system),
      identifier_(system),
      hist_(),
      speciesId_(-1),
      mutatorPtr_(0),
      outputFileName_(),
      initialWeightFileName_("0"),
      stepCount_(0)
   {  setClassName("SingleMicelleUmbrellaSamplingMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void SingleMicelleUmbrellaSamplingMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "atomType", atomTypeId_);
      // Cast the Species to LinearSG
      speciesPtr_ = dynamic_cast<LinearSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be LinearSG");
      }
      capacity_ = speciesPtr_->capacity()+1; 
      mutatorPtr_ = &speciesPtr_->mutator();
      weights_.allocate(capacity_);
      oldIdentities_.allocate(capacity_);
      read<int>(in, "aggregateCutoff", aggregateCutoff_);
      read<int>(in, "flips", flipsPerMove_);
      read<double>(in, "cutoff", cutoff_);
      read<int>(in, "histMin", histMin_);
      read<int>(in, "histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      read<int>(in, "upperLimit", uLimit_);
      read<int>(in, "lowerLimit", lLimit_);
      read<std::string>(in, "clusterHist",clusterHistName_);
      read<std::string>(in, "outputFileName",outputFileName_);
      readOptional<std::string>(in, "initialWeights",initialWeightFileName_);
      std::ifstream weightFile;
      if (initialWeightFileName_!="0") {
        system().fileMaster().open(initialWeightFileName_, weightFile);
        int n;
        double m;
        while (weightFile >> n >>m) {
          weights_[n]= m;
        }
      } else {
        for (int x = 0; x < capacity_; ++x) {
          weights_[x]=0;
        }
      }
      system().fileMaster().openOutputFile(clusterHistName_, clusterFile_);
      flipsAttempted_=0;
      flipsAccepted_=0;
      MicelleRejections_=0;
   }
   /*
   * Load state from an archive.
   */
   void SingleMicelleUmbrellaSamplingMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "atomType", atomTypeId_);
      // Cast the Species to LinearSG
      speciesPtr_ = dynamic_cast<LinearSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a LinearSG");
      }
      loadParameter<int>(ar, "aggregateCutoff", aggregateCutoff_);
      loadParameter<int>(ar, "flips", flipsPerMove_);
      loadParameter<double>(ar, "cutoff", cutoff_);
      loadParameter<int>(ar, "histMin", histMin_);
      loadParameter<int>(ar, "histMax", histMax_);
      loadParameter<int>(ar, "UpperLimit", uLimit_);
      loadParameter<int>(ar, "LowerLimit", lLimit_);
      mutatorPtr_ = &speciesPtr_->mutator();
      ar & weights_;
   }
   

   /*
   * Save state to an archive.
   */
   void SingleMicelleUmbrellaSamplingMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & aggregateCutoff_;
      ar & flipsPerMove_;
      ar & cutoff_;
      ar & histMin_;
      ar & histMax_;
      ar & uLimit_;
      ar & lLimit_;
      ar & weights_;
   }
     

   Molecule& SingleMicelleUmbrellaSamplingMove::randomSGMolecule(int speciesId, int nSubType, int flipType, bool atBoundary)
   {  // If the system is not on the bound any random molecule will do
      if (!atBoundary) {
        return system().randomMolecule(speciesId_);
      }
      // Otherwise choose a random molecule of the necessary subtype
      int moleculeId,nMol,index,type;
      int count = 0;
      nMol = system().nMolecule(speciesId);
      if (nMol <= 0) {
        Log::file() << "Number of molecules in species " << speciesId
                     << " = " << nMol << std::endl;
        UTIL_THROW("Number of molecules in species <= 0");
      }
      index = simulation().random().uniformInt(0, nSubType);
      for (int i=0; i<nMol; ++i) {
        type = speciesPtr_->mutator().moleculeStateId(system().molecule(speciesId, i));
        if (type==flipType) {
          if (count==index) {
            moleculeId = i;
            return system().molecule(speciesId, moleculeId);
          }
          count += 1;
        }
      }
      // If the code loops through all molecules and can't find the indexed molecule
      UTIL_THROW("No molecule selected");
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool SingleMicelleUmbrellaSamplingMove::move() 
   { 
      incrementNAttempt();
      bool accept;
      System::MoleculeIterator molIter;
      system().begin(speciesId_, molIter);
      int moleculeId;
      Cluster thisCluster;
      int clusterSize, nClusters;
      double comboPrefactor = 1.0;
      int flipType = -1.0;
      int flipTypeCapacity;
      bool atBoundary;
      bool wasAtUpperBound;
      bool wasAtLowerBound;
      // Test if the system starts with one cluster
      identifier_.identifyClusters();
      nClusters=identifier_.nCluster();
      int clusterCount = 0; 
      for (int i = 0; i < nClusters; i++) {
        thisCluster=identifier_.cluster(i);
        clusterSize=thisCluster.size();
        if (clusterSize > aggregateCutoff_) {
          clusterCount=clusterCount+1;
        }
      }
      if (clusterCount>2) {
         UTIL_THROW("More than one cluster at move start!!!");
      }
      // Record chemical identity of each molecule before the flip chain begins
      for ( int i =0; i < capacity_-1; ++i) {
        Molecule& currentMolecule = system().molecule(speciesId_,i);
        oldIdentities_[i] = speciesPtr_->mutator().moleculeStateId(currentMolecule);
      }

      // Perform the desired number of flipping moves
      for (int i =0; i < flipsPerMove_; ++i) {
        flipsAttempted_ +=1;
        comboPrefactor = 1.0;
        atBoundary=false;
        wasAtLowerBound=false;
        wasAtUpperBound=false;
        int oldStateCount = mutatorPtr_->stateOccupancy(0);
        int oldState = mutatorPtr_->stateOccupancy(0);
        // The next two if statements allow the code to always perform a flip away from the boundary when at the upper 
        // or lower limit for number of allowed copolymers. This is accomplished by always flipping a molecule of the given type,
        // but correcting for this via a combinatorial prefactor
        if  (oldStateCount == uLimit_) {
          atBoundary = true;
          flipType = 0;
          flipTypeCapacity =  mutatorPtr_->stateOccupancy(0);
          wasAtUpperBound=true;
        }
        if (oldStateCount == lLimit_) {
          atBoundary = true;
          wasAtLowerBound=true;
          flipType = 1;
          flipTypeCapacity =  mutatorPtr_->stateOccupancy(1); 
        }
        // Gets a random molecule of the needed type
        Molecule& molecule = randomSGMolecule(speciesId_, flipTypeCapacity, flipType, atBoundary);
        #ifndef INTER_NOPAIR
        // Calculate pair energy for the chosen molecule
        double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
        #endif
        // Toggle state of the molecule
        int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
        int newStateId = (oldStateId == 0) ? 1 : 0;
        speciesPtr_->mutator().setMoleculeState(molecule, newStateId);
        #ifdef INTER_NOPAIR 
        bool   accept = true;
        #else //ifndef INTER_NOPAIR
        // Recalculate pair energy for the molecule
        double newEnergy = system().pairPotential().moleculeEnergy(molecule);
        // Calculate the probability of acceptance and determine if the move is accepted or rejected
        int    newState = mutatorPtr_->stateOccupancy(0);
        // If the system was at either the upper or lower bound calculate the combinatorial prefactor necessary for detailed balance
        if (wasAtLowerBound) {
          comboPrefactor = (double)(capacity_-1-lLimit_)/(double)(capacity_-1);
        }
        if (wasAtUpperBound) {
          comboPrefactor = (double)(uLimit_)/(double)(capacity_-1);
        }
        // Calculate changees in 
        double oldWeight = weights_[oldState];
        double newWeight = weights_[newState];
        double oldWeightSG = speciesPtr_->mutator().stateWeight(oldStateId);
        double newWeightSG = speciesPtr_->mutator().stateWeight(newStateId);
        double ratio  = boltzmann(newEnergy - oldEnergy)*boltzmann(newWeight-oldWeight)*newWeightSG/oldWeightSG*comboPrefactor;
           accept = random().metropolis(ratio);
        #endif
        if (accept) {
        
        flipsAccepted_ +=1;
        } else {
        // Revert chosen molecule to original state
          speciesPtr_->mutator().setMoleculeState(molecule, oldStateId);
        }
      }

      // Test if a second cluster has formed
      identifier_.identifyClusters();
      nClusters=identifier_.nCluster();
      clusterCount = 0; 
      for (int i = 0; i < nClusters; i++) {
        thisCluster=identifier_.cluster(i);
        clusterSize=thisCluster.size();
        if (clusterSize > aggregateCutoff_) {
          clusterCount=clusterCount+1;
        }
      }
      //If so reject the sequence, record the failed move, and return back to the prior state 
      if (clusterCount > 1) {
        MicelleRejections_ +=1;
        accept = false;
        system().begin(speciesId_, molIter);
        for (int i ; i<capacity_-1; ++i) {
          Molecule currentMolecule = system().molecule(speciesId_,i);
          speciesPtr_->mutator().setMoleculeState(currentMolecule, oldIdentities_[i]);
        }
        hist_.clear();
        identifier_.identifyClusters();
        for (int i = 0; i < identifier_.nCluster(); i++) {
          hist_.sample(identifier_.cluster(i).size());
        }
        int min = hist_.min();
        int nBin = hist_.nBin();
        for (int i = 0; i < nBin; ++i) {
           clusterFile_ << Int(i + min) << "  "
                   <<  Dbl(double(hist_.data()[i])) << "\n";
        }
      } else {
        accept = true;
        incrementNAccept();
        hist_.clear();
        for (int i = 0; i < identifier_.nCluster(); i++) {
          hist_.sample(identifier_.cluster(i).size());
        }
        int min = hist_.min();
        int nBin = hist_.nBin();
        for (int i = 0; i < nBin; ++i) {
           clusterFile_ << Int(i + min) << "  "
                   <<  Dbl(double(hist_.data()[i])) << "\n";
        }
      }       

      return accept;

   }
   // Brief output of all the flips and how many were rejectred for creating 2 or more micelles
   void SingleMicelleUmbrellaSamplingMove::output()
   {    clusterFile_.close();
        std::ofstream outputFile;
        system().fileMaster().openOutputFile("clusterAcceptanceStats",outputFile);
        outputFile << flipsAttempted_ <<"  " << flipsAccepted_ << "  " << MicelleRejections_;
        outputFile.close();
        std::string fileName = outputFileName_; 
        fileName += ".dat";
        system().fileMaster().openOutputFile(fileName, outputFile);
        for (int i = 0; i < capacity_; i++) {
           outputFile << i << "   " <<  weights_[i]<<std::endl;
        }
        outputFile.close();
   }
   

}
