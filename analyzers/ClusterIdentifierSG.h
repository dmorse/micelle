#ifndef MCMD_CLUSTER_IDENTIFIER_SG_H
#define MCMD_CLUSTER_IDENTIFIER_SG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/Cluster.h>       // member template argument
#include <mcMd/analyzers/system/ClusterLink.h>   // member template argument
#include <mcMd/neighbor/CellList.h>              // member
#include <simp/boundary/Boundary.h>              // argument (typedef)
#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template
#include <util/containers/GStack.h>              // member template

namespace Simp {
   class Species;
}

namespace McMd
{

   class System;

   using namespace Util;
   using namespace Simp;

   /**
   * Identifies clusters of molecules, such as micelles.
   */
   class ClusterIdentifierSG
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClusterIdentifierSG(System &system);
   
      /**
      * Destructor.
      */
      ~ClusterIdentifierSG();
   
      /** 
      * Initialize state and allocate required memory (call once).
      *
      * \param speciesId index of species in clusters
      * \param atomTypeId typeId of atoms in micelle core
      * \param cutoff pair distance cutoff
      * \param subType index of SG molecule subtype to include in clusters
      */
      virtual 
      void initialize(int speciesId, int atomTypeId, double cutoff, int subType);
   
      /**
      * Identify all clusters (main operation).
      */
      void identifyClusters();

      /**
      * Get number of clusters.
      */ 
      int nCluster() const
      {  return clusters_.size(); }

      /**
      * Get a specific ClusterLink, by id of the associated molecule.
      *
      * \param moleculeId molecule index.
      */ 
      ClusterLink& link(int moleculeId)
      {  return links_[moleculeId]; }

      /**
      * Get a specific cluster, indexed in the order identified.
      *
      * The id argument of this function is a consecutive array
      * index, with 0 <= id < nCluster, which need not be equal 
      * to the cluster identifier returned by Cluster::id().
      *
      * \param id cluster array index, 0 <- id < nCluster.
      */ 
      Cluster& cluster(int id)
      {  return clusters_[id]; }

      /**
      * Return true if valid, or throw Exception otherwise.
      */
      bool isValid() const;

   private:

      /// Array of ClusterLink objects, indexed by molecule id.
      DArray<ClusterLink> links_;

      /// Growable array of clusters.
      GArray<Cluster> clusters_;

      /// Work stack of unprocessed cluster links.
      GStack<ClusterLink>  workStack_;

      /// CellList of atoms of the specified species and atom type.
      CellList cellList_;

      /// Pointer to parent System.
      System* systemPtr_;

      /// Molecule species type id
      int speciesId_;

      /// Atom typeId of core atoms
      int atomTypeId_;

      /// Cutoff distance for touching cores
      double cutoff_;

      /// Subtype index for molecules that can form a cluster. 
      int subTypeId_;

      // Private functions

      /**
      * Return parent system by reference.
      */
      System& system()
      {  return *systemPtr_; }

      /**
      * Pop and process top ClusterLink in the workStack.
      */
      void processNextMolecule(Cluster& cluster);

   };

}
#endif
