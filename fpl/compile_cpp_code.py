import fenics as fe
import os

srcdir = str(os.path.dirname(os.path.realpath(__file__))+"/src/")
header_file = open(srcdir+"/Particles.h", "r")
code = header_file.read()
header_file.close()

compiled_module = fe.compile_extension_module(
    code=code, source_directory=srcdir,
    sources=["Particles.cpp"],
    include_dirs=[".",os.path.abspath(srcdir)],
    additional_declarations="""
%feature("notabstract") Particles;
%template(SizetVector) std::vector<std::size_t>;
%template(DoubleVector) std::vector<double>;
"""
)

CalcF = compiled_module.Particles.CalcF

