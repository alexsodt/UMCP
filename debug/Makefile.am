AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov


bin_PROGRAMS = hd.dbg makeHexSquare.dbg min.dbg subdivide.dbg scale.dbg icosahedron.dbg makeBud.dbg makeBump.dbg makeHexCylinder.dbg writeBiasPotential.dbg writeCylRTraj.dbg rd.dbg makePore.dbg join.dbg

AM_CPPFLAGS = -include ../config.h -include ../src/local_config.h $(AM_INCLUDE) 
AM_CFLAGS = -g
AM_CXXFLAGS = -g

vpath %.C ../src
vpath %.c ../src
vpath %.cpp ../src
vpath %.h ../src


complex_SRCS = NBAR.C dimer.C crowder.C
util_SRCS    = mutil.C util.C

interp_SRCS = interpLib.C grad.C pboxing.C irr_kernel.C uv_map.C gauss.C modes.C write_structures.C min_distance.cpp GJK.cpp gjk_gudrun.cpp GJK_min2.cpp functions.cpp NearestSimplexR2.cpp NearestSimplexR.cpp  subdivide3.C fast_mm.C random_global.C height.C pcomplex.C minimize.C gsl_random_globals.C dynamics.C restart.C p_p.C sparse.C global_boxing.C lipid_density.C tachyon.C npt.C lipid_composition.C vpass.C get_regions.C parallel.C input.C l-bfgs.C srd.C alignSet.C compilation_time.C meshCollisionLib.C sans.C spline.C globals.C createAllAtom.C dcd.C pdb.C simulation.C

rd_SRCS = fpr_subroutines/Faddeeva.cpp fpr_subroutines/passocF.cpp fpr_subroutines/pirr_pfree_ratio_psF.cpp 

hull_SRCS = dovoronoi.c clusterv.C

rd_dbg_SOURCES        = rd.C $(rd_SRCS) $(interp_SRCS) $(complex_SRCS) $(util_SRCS)
hd_dbg_SOURCES        = hd.C $(interp_SRCS) $(complex_SRCS) $(util_SRCS)
min_dbg_SOURCES       = min.C $(interp_SRCS) $(complex_SRCS) $(util_SRCS)
subdivide_dbg_SOURCES = subdivide.C $(interp_SRCS) $(complex_SRCS) $(util_SRCS)
scale_dbg_SOURCES	  = scale.C  $(interp_SRCS) $(complex_SRCS) $(util_SRCS)
icosahedron_dbg_SOURCES = icosahedron.C $(util_SRCS)
join_dbg_SOURCES	= join.C $(interp_SRCS) $(complex_SRCS) $(util_SRCS) $(hull_SRCS) getCylinderMesh.C
makeBud_dbg_SOURCES   = makeBud.C $(hull_SRCS) $(util_SRCS)
makePore_dbg_SOURCES   = makePore.C $(hull_SRCS) $(util_SRCS) getCylinderMesh.C
makeBump_dbg_SOURCES  = makeBump.C $(hull_SRCS) $(util_SRCS)
makeHexCylinder_dbg_SOURCES = makeHexCylinder.C $(hull_SRCS) $(util_SRCS)
makeHexSquare_dbg_SOURCES = makeHexSquare.C $(hull_SRCS) $(util_SRCS)
writeBiasPotential_dbg_SOURCES = writeBiasPotential.C $(util_SRCS)
writeCylRTraj_dbg_SOURCES      = writeCylRTraj.C $(util_SRCS) 


rd_dbg_LDADD		  = $(AM_LDADD) -lstdc++ 
hd_dbg_LDADD              = $(AM_LDADD) -lstdc++
join_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) -lstdc++
min_dbg_LDADD             = $(AM_LDADD) -lstdc++
subdivide_dbg_LDADD       = $(AM_LDADD) -lstdc++
scale_dbg_LDADD           = $(AM_LDADD) -lstdc++
icosahedron_dbg_LDADD     = $(AM_LDADD) -lstdc++
makeBud_dbg_LDADD         = $(AM_LDADD) $(AM_QHULL) -lstdc++
makePore_dbg_LDADD         = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeBump_dbg_LDADD        = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexCylinder_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexSquare_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeBiasPotential_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeCylRTraj_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++

compilation_time.o: .FORCE
.PHONY: .FORCE

