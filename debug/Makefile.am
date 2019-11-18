AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = hd.dbg makeHexSquare.dbg min.dbg subdivide.dbg scale.dbg icosahedron.dbg makeBud.dbg makeBump.dbg makeHexCylinder.dbg writeBiasPotential.dbg writeCylRTraj.dbg makePore.dbg join.dbg frc_sandbox.dbg shiftMesh.dbg

AM_CPPFLAGS = -include ../config.h -include ../src/local_config.h $(AM_INCLUDE) 
AM_CFLAGS = -g -pg
AM_CXXFLAGS = -g -pg

SRC_DIR = ../src

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)
vpath %.h  $(SRC_DIR)
vpath %.cu $(SRC_DIR)

.cu.o:
	$(NVCC) -o $*.o -c $(SRC_DIR)/$*.cu

complex_SRCS = NBAR.C dimer.C crowder.C
util_SRCS    = mutil.C util.C
amoeba_SRCS  = amoeba_code/amotry.C amoeba_code/amoeba.C amoeba_code/nrutil.C 
rd_SRCS = fpr_subroutines/Faddeeva.cpp fpr_subroutines/passocF.cpp fpr_subroutines/pirr_pfree_ratio_psF.cpp fpr_subroutines/2D.cpp fpr_subroutines/DDmatrixcreate.cpp fpr_subroutines/DDpsur.cpp fpr_subroutines/TBLnorm.cpp fpr_subroutines/TBLpirr.cpp fpr_subroutines/TBLsur.cpp fpr_subroutines/fpir.cpp fpr_subroutines/fsur.cpp fpr_subroutines/integrator.cpp fpr_subroutines/fnorm.cpp fpr_subroutines/sizelookup.cpp rd.C fpr_subroutines/DDpirr_pfree_ratio_ps.cpp fpr_subroutines/pirr.cpp fpr_subroutines/pnorm.cpp

interp_SRCS = interpLib.C grad.C pboxing.C irr_kernel.C uv_map.C gauss.C modes.C write_structures.C min_distance.cpp GJK.cpp gjk_gudrun.cpp GJK_min2.cpp functions.cpp NearestSimplexR2.cpp NearestSimplexR.cpp  subdivide3.C fast_mm.C random_global.C height.C pcomplex.C minimize.C gsl_random_globals.C dynamics.C restart.C p_p.C sparse.C global_boxing.C lipid_density.C tachyon.C npt.C lipid_composition.C vpass.C get_regions.C parallel.C input.C l-bfgs.C srd.C alignSet.C compilation_time.C meshCollisionLib.C sans.C spline.C globals.C createAllAtom.C dcd.C pdb.C simulation.C init.C 3d_interp.C fitRho.C M_matrix.C gather.C io_mol_read.C $(amoeba_SRCS) $(rd_SRCS)


hull_SRCS = dovoronoi.c clusterv.C

if CUDA
cuda_SRCS = cuda.cu
else
cuda_SRCS =
endif

hd_dbg_SOURCES        = hd.C $(interp_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
frc_sandbox_dbg_SOURCES        = frc_sandbox.C $(interp_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
shiftMesh_dbg_SOURCES        = shiftMesh.C $(interp_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
min_dbg_SOURCES       = min.C $(interp_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) 
subdivide_dbg_SOURCES = subdivide.C $(interp_SRCS)  $(hull_SRCS) $(complex_SRCS) $(util_SRCS)
scale_dbg_SOURCES	  = scale.C  $(interp_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS)
icosahedron_dbg_SOURCES = icosahedron.C $(util_SRCS)
join_dbg_SOURCES	= join.C $(interp_SRCS) $(complex_SRCS) $(util_SRCS) $(hull_SRCS) getCylinderMesh.C
makeBud_dbg_SOURCES   = makeBud.C $(hull_SRCS) $(util_SRCS)
makePore_dbg_SOURCES   = makePore.C $(hull_SRCS) $(util_SRCS) getCylinderMesh.C
makeBump_dbg_SOURCES  = makeBump.C $(hull_SRCS) $(util_SRCS)
makeHexCylinder_dbg_SOURCES = makeHexCylinder.C $(hull_SRCS) $(util_SRCS)
makeHexSquare_dbg_SOURCES = makeHexSquare.C $(hull_SRCS) $(util_SRCS)
writeBiasPotential_dbg_SOURCES = writeBiasPotential.C $(util_SRCS)
writeCylRTraj_dbg_SOURCES      = writeCylRTraj.C $(util_SRCS) 


hd_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
frc_sandbox_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
shiftMesh_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
join_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) -lstdc++
min_dbg_LDADD             = $(AM_LDADD) $(AM_QHULL) -lstdc++
subdivide_dbg_LDADD       = $(AM_LDADD) $(AM_QHULL) -lstdc++
scale_dbg_LDADD           = $(AM_LDADD) $(AM_QHULL) -lstdc++
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

