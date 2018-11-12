using Pkg
Pkg.activate(".")
using MCMphotolysis

jvals = j_oldpars("testM4", output = false)
jvals = j_oldpars("testM4")
