../../optimized/icosahedron.opt > icos.mesh
../../optimized/subdivide.opt icos.mesh
../../optimized/subdivide.opt subdiv.mesh
../../optimized/min.opt subdiv.mesh
../../optimized/scale.opt min.mesh 300
mv scale.mesh sphere.mesh
#mv sphere.inp_force sphere.inp

