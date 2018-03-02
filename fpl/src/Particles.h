#ifndef __FPL_PARTICLES_H
#define __FPL_PARTICLES_H

#include <string>
#include <utility>
#include <vector>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/PETScVector.h>

#include <dolfin/log/log.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/geometry/BoundingBoxTree.h>

namespace dolfin
{
  /// This class represents a collection of discrete particles.
  /// It can be used with assembly routines to couple to Functions
  class Particles
  {
  public:
    Particles(int NP, int dim);
    static void CalcF(GenericVector * px, GenericVector * pv, GenericVector * pr,
		      std::size_t nploc,

		      Function &fluid_u,
		      GenericVector * f_f2p, GenericVector * f_p2f);
  };

}

#endif
