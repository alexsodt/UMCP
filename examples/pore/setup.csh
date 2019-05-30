# Syntax: makeBud NX LX neckRadius neckHeight SphereRadius
makePore.opt 6 180 50 60 30
subdivide.opt fusion.mesh
min.opt subdiv.mesh
min.opt min.mesh
min.opt min.mesh
min.opt min.mesh
mv min.mesh pore.mesh
