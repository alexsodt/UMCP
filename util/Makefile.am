AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = getPoreInfo.dbg 

AM_CPPFLAGS = -include ../config.h -include ../src/local_config.h $(AM_INCLUDE) 
AM_CFLAGS = -g
AM_CXXFLAGS = -g

SRC_DIR = ../src

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)
vpath %.h  $(SRC_DIR)
vpath %.cu $(SRC_DIR)


getPoreInfo_dbg_SOURCES = getPoreInfo.C util.C

getPoreInfo_dbg_LDFLAGS	  = $(AM_LDADD) -lstdc++

