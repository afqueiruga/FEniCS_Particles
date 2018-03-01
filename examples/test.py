import fenics as fe
import mshr

mesh = RectangleMesh(-2,0,2,1 ,40,40, "right/left")
boundaries = FacetFunction("uint", mesh)
boundaries.set_all(0)
