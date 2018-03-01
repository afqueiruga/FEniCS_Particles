
#include "Particles.h"
#include <dolfin/fem/fem_utils.h>
#include <dolfin/geometry/Point.h>

using namespace dolfin;

double norm(double * x, int dim) {
  double ac = 0.0;
  int i;
  for(i=0;i<dim;i++) ac+=x[i]*x[i];
  return sqrt(ac);
}
void Kernel_F2P(int dim,
		double * px, double *pv, double *pr,
		double *u, double *f) {
  double vr[dim];
  for(int i=0;i<dim;i++) {
    vr[i]=pv[i]-u[1+i]/u[0]; // v=g/rho
  }
  double PIR2 = M_PI*pr[0]*pr[0];
  double K = 0.5*u[0]*0.47*PIR2;
  for(int i=0;i<dim;i++) {
    f[i]= -K*norm(vr,dim)*vr[i];
  }
}

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
    //Get the local particle values
    for(int j=0;j<dim;j++) v_dof_idx[j] = dim*i+j;
    s_dof_idx[0] = i;
    px->get_local(p_x,dim,v_dof_idx);
    pv->get_local(p_v,dim,v_dof_idx);
    pr->get_local(p_r,1,s_dof_idx);
    // Get the cell that the particle is in
    // TODO: Optimize with a table. The particles need a cell table
    // anyways for distribution.
    Point p_Px(dim,p_x);
    unsigned int cix = f_bbox->compute_first_entity_collision(p_Px);
    // Skip to the next particle, this probably means it's not in the mesh anymore. Oh well.
    if(cix >= f_U->mesh()->num_cells()) continue;
    // Grab the data for the local cell
    const Cell f_cell(*f_U->mesh(),cix);
    f_cell.get_vertex_coordinates(f_vertex_coordinates);
    // dolfin::ArrayView<const int>
    // const std::vector<dolfin::la_index> f_dofs = f_dmap->cell_dofs(cix);
    const int c_orient = f_cell.orientation();
    // Eval fluid at particle position
    fluid_u(f_u,p_x[0],p_x[1]); //TODO: Optimize call with cell info

    double kernel_force[dim];
    Kernel_F2P(dim, p_x,p_v,p_r, f_u.data(), kernel_force);
    
    if(f_f2p) {
      f_f2p->add_local(kernel_force,dim, v_dof_idx);
    }
  }
}
