#ifndef MCMD_LINEAR_SG_H
#define MCMD_LINEAR_SG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesMutator.h"         // base class
#include <simp/species/Linear.h>    // base class
#include <util/containers/DArray.h> // member template
#include <util/containers/Pair.h>   // member

namespace McMd
{

   using namespace Util;

   /**
   * A mutable linear polymer, for semigrand ensemble.
   *
   * A LinearSG molecule is a mutable Linear chain that can be in
   * either of two states, and in which atom types can be toggled
   * between two corresponding sequences of values. 
   *
   * \ingroup McMd_Species_Module
   */
   class LinearSG : public Linear, public SpeciesMutator
   {
   
   public:
  
      /**
      * Constructor.
      */
      LinearSG();
   
      /** 
      * Destructor.
      */
      virtual ~LinearSG();
  
      /**
      * Set the type of all atoms in the molecule.
      *
      * \param molecule reference to molecule of interest.
      * \param stateId  atom type id.
      */
      virtual void setMoleculeState(Molecule& molecule, int stateId);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);


   protected:
   
      /**
      * Read nAtom, a pair of atom type ids and weightRatio.
      *
      * \param in input stream
      */
      virtual void readSpeciesParam(std::istream &in);
   
      /**
      * Load species structure from an Archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

      /**
      * Return the same type for any particle in any chain.
      *
      * \param index atom index, in range 0,...,nAtom_ - 1
      * \return atom type index
      */
      virtual int calculateAtomTypeId(int index) const;
   
      /**
      * Return same bond type for any bond in any chain.
      *
      * \param index bond index, in range 0,...,nAtom_ - 1
      * \return bond type index
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef INTER_ANGLE
      /**
      * Return same angle type for any angle in any chain.
      *
      * \param index  angle index, in range 0, ..., nAngle_ - 1
      * \return       angle type index
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Return same dihedral type for any dihedral in any chain.
      *
      * \param index  dihedral index, in range 0, ..., nDihedral_ - 1
      * \return       dihedral type index
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif


   private:

      // A Pair of atom type indexes.
      // Atoms in a molecule with stateId = i have type typeIds_[i]
      Pair<int> typeIds_;

      /// Ratio of statistical weights for typeId[0]/typeId[1]
      double weightRatio_;
      
      // Used for setting the molecule state to type 1 or type 2
      DArray<int> beadIdentities_;

      // Array which has the atomtype for each bead in species 1
      DArray<int> beadTypeIds1_;

      // Array which has the atomtype for each bead in species 2
      DArray<int> beadTypeIds2_; 

      // Maximum placement attempts when attempting genergate molecules
      // May no longer be needed
      static const int maxPlacementAttempts_ = 500;

   };
   
} 
#endif
