#include "Particles.h"
using namespace dolfin;

Particles::Particles(int NP, int dim) {

}

void Particles::CalcF(GenericVector * px, GenericVector * pv, GenericVector * pr,
		      std::size_t nploc,
		      Function &fluid_u,
		      GenericVector * f_f2p, GenericVector * f_p2f)
{

  
  std::shared_ptr<const FunctionSpace> f_U = fluid_u.function_space();
  std::shared_ptr<BoundingBoxTree> f_bbox = f_U->mesh()->bounding_box_tree();
  std::shared_ptr<const GenericDofMap> f_dmap = f_U->dofmap();

  // Extract the subspace
  std::shared_ptr<const FiniteElement> f_elem = f_U->element();
  std::vector<std::size_t> f_g_ix; f_g_ix.push_back(1);
  std::shared_ptr<const FunctionSpace> f_g_U = f_U->extract_sub_space(f_g_ix);
  std::shared_ptr<const FiniteElement> f_g_elem = f_elem->extract_sub_element(f_g_ix);

  
  // Temporary
  int dim = 2;
  double p_x[dim], p_v[dim], p_r[1];
  Array<double> f_u(2+dim);
  std::vector<double> f_vertex_coordinates;
  la_index v_dof_idx[dim];
  la_index s_dof_idx[1];
  for(std::size_t i=0;i<nploc;i++) {
    //Get the local field values

    for(int j=0;j<dim;j++) v_dof_idx[j] = dim*i+j;
    s_dof_idx[0] = i;
    px->get_local(p_x,dim,v_dof_idx);
    pv->get_local(p_v,dim,v_dof_idx);
    pr->get_local(p_r,1,s_dof_idx);
  }
}
