import fenics as fe
from compile_cpp_code import CalcF
import numpy as np

class Particles():
    def __init__(self, Np, dim, r=1.0):
        self.Np = Np
        self.dim = dim
        self.x = fe.Vector()
        self.x.init(Np*dim)
        self.x.zero()
        self.v = fe.Vector()
        self.v.init(Np*dim)
        self.v.zero()
        self.r = fe.Vector()
        self.r.init(Np)
        self.r[:] = r
        
    def init_grid(self, bx,by,tx,ty, dev = 0.2):
        Nside = int(self.Np ** (1.0/self.dim))
        lx = self.x.get_local()
        for i in xrange(self.Np):
            lx[self.dim*i:self.dim*(i+1)] = [
                (tx-bx)/Nside*np.floor(i/Nside)+bx,
                (ty-by)/Nside*np.remainder(i,Nside)+by] + \
                np.random.uniform(-dev,dev,self.dim)
        self.x.set_local(lx)
        self.x.apply("insert")

    def write_vtk(self, fname):
        pass
    
    def euler_step(self, u, DT):
        "Take a step using forward euler. Good enough for tracers."
        F = fe.Vector()
        F.init(self.Np*self.dim)
        F.zero()
        CalcF(self.x,self.v,self.r,
              self.Np,
              u,
              F,None)
        self.x += DT*self.v
        self.v += DT*F
        
    def RK_Field(self):
        "Returns an RK_Field object for use with afqsrungekutta"
        return None
