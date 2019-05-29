# Syntax: makeBud NX LX neckRadius neckHeight SphereRadius
makePore.opt 6 180 50 80 30
subdivide.opt fusion.mesh
min.opt subdiv.mesh
mv min.mesh fusion.mesh
