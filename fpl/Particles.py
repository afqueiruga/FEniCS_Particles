import fenics as fe
from compile_cpp_code import CalcF

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
