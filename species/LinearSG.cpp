/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/neighbor/CellList.h>
#include <mcMd/chemistry/Molecule.h>
#include "LinearSG.h"
#ifdef UTIL_MPI
#include <mcMd/simulation/McMd_mpi.h>
#endif

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LinearSG::LinearSG()
    : Linear(),
      SpeciesMutator(),
      beadTypeIds1_(),
      beadTypeIds2_()
   {
      setMutatorPtr(this);
   }

   /*
   * Destructor.
   */
   LinearSG::~LinearSG()
   {}

   /*
   * Read atom structure and two sets of atom type ids.
   */
   void LinearSG::readSpeciesParam(std::istream& in)
   {
      read<int>(in,"nAtom", nAtom_);
      nBond_ = nAtom_ - 1;
      #ifdef INTER_ANGLE
      nAngle_ = nAtom_ - 2;
      #endif
      #ifdef INTER_DIHEDRAL
      if (nAtom_ > 3)
         nDihedral_ = nAtom_ - 3;
      else
         nDihedral_ = 0;
      #endif
      buildLinear();

      read<Pair <int> >(in, "typeIds", typeIds_);
      beadTypeIds1_.allocate(nAtom_);
      beadTypeIds2_.allocate(nAtom_);
      readDArray<int>(in, "identities1", beadTypeIds1_, nAtom_);
      readDArray<int>(in, "identities2", beadTypeIds2_, nAtom_);

      read<double>(in, "weightRatio", weightRatio_);

      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   void LinearSG::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar,"nAtom", nAtom_);
      nBond_  = nAtom_ - 1;
      #ifdef SIMP_ANGLE
      hasAngles_ = 0;
      loadParameter<int>(ar,"hasAngles", hasAngles_, false);
      if (hasAngles_) {
         nAngle_ = nBond_ - 1;
         if (nAngle_ > 0) {
            loadParameter<int>(ar,"angleType", angleType_);
         }
      } else {
         nAngle_ = 0;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      hasDihedrals_ = 0;
      loadParameter<int>(ar,"hasDihedrals", hasDihedrals_, false);
      if (hasDihedrals_) {
         if (nAtom_ > 3) {
            nDihedral_ = nAtom_ - 3;
         } else {
            nDihedral_ = 0;
         }
         if (nDihedral_ > 0) {
            loadParameter<int>(ar, "dihedralType", dihedralType_);
         }
      } else {
         nDihedral_ = 0;
      }
      #endif
      buildLinear();
      loadParameter<Pair <int> >(ar, "typeIds", typeIds_);
      
      ar & beadTypeIds1_;
      ar & beadTypeIds2_;

      loadParameter<double>(ar, "weightRatio", weightRatio_);
      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   /*
   * Save internal state to an archive.
   */
   void LinearSG::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar << nAtom_;
      ar << typeIds_;
      ar << beadTypeIds1_;
      ar << beadTypeIds2_;
      ar << weightRatio_;
      #ifdef SIMP_ANGLE
      Parameter::saveOptional(ar, hasAngles_, hasAngles_);
      if (hasAngles_ && nAngle_ > 0) {
         ar << angleType_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      Parameter::saveOptional(ar, hasDihedrals_, hasDihedrals_);
      if (hasDihedrals_ && nDihedral_ > 0) {
         ar << dihedralType_;
      }
      #endif
   }

   /*
   * Return NullIndex for every atom.
   * Set initial typeIds
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAtomTypeId(int index) const
   { return 0; }

   /*
   * Return 0 for every bond.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateBondTypeId(int index) const
   { return 0; }

   #ifdef INTER_ANGLE
   /*
   * Return 0 for every angle.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAngleTypeId(int index) const
   { return 0; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Return 0 for every dihedral.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateDihedralTypeId(int index) const
   { return 0; }
   #endif

   /*
   * Change the type of a specific molecule.
   */
   void LinearSG::setMoleculeState(Molecule& molecule, int stateId)
   {
      int nAtom  = molecule.nAtom();
      if (stateId == 0) {
        beadIdentities_ = beadTypeIds1_;
      } else {
        beadIdentities_ = beadTypeIds2_;
      }
      for (int i = 0; i < nAtom; ++i) {
        molecule.atom(i).setTypeId(beadIdentities_[i]);
      }
      setMoleculeStateId(molecule, stateId);
   }
}
