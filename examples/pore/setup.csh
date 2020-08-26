# approximately matches boundary conditions.. sufficient density (6)
# Lx Lz Ly
../../optimized/makeHexSquare.opt 6 190 195 190

# mesh pt 72 is near the middle. I am going to change this so that it doesn't matter..
set pt1 = 72
set pt2 = 72
set rad = 45
set H = 100 

../../debug/join.dbg planar.mesh $pt1 planar.mesh $pt2 $rad $H auto
../../optimized/hd.opt min.inp mesh=join.mesh jobname=join1 nmin=2 movie=yes outputMesh=yes > hd_opt.out 
mv join1.mesh pore.mesh
