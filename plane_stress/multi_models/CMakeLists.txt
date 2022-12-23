## CMAKE TO CREATE A SHARED LIBRARY FOR ABQ USERSUBROUTINES
## Include EIGEN C++ linearalgebra library (header based)

###########################################################
# DISROOT - Where to save the library file (same adress in abq_v6.env)
# ABQROOT - aba_to_c.h location
# EIGENROOT - Location where EIGEN folder is
# ---------------------------------------------------------
#  add_library(Library_name STATIC subroutine.cpp support.cpp ...)
###########################################################

cmake_minimum_required(VERSION 2.6)
project (UserSubLib)

## set(BLASROOT "/usr/local/lib/libblas.a")
set(EIGENROOT "/home/vasudevan/PhD/Eigen")
set(ABQROOT "/usr/local/abaqus/6.14-2/6.14-2/code/include")
set(DISROOT "/home/vasudevan/PhD/Code/UMAT/plane_stress_multi_models/new_ver")

#main.cpp

add_library(UserSubLib STATIC Polynomial.cpp Hencky.cpp dashpot.cpp elastic_branch.cpp) # MODIFY SOURCE FILES HERE

### set (CMAKE_CXX_FLAGS "-ffast-math -libblas -fopenmp -fPIC -DPIC -w -Wuninitialized \
###	    -DHKS_LINUX -DLINUX -DHAS_BOOL -DFOR_TRAIL -DHKS_DEBUG \
###	    -D_BSD_TYPES -D_BSD_SOURCE -D_GNU_SOURCE \
###	    -D_XOPEN_SOURCE_EXTENDED -DEMULATE_EXCEPTIONS=0 \
###	    -DGL_GLEXT_PROTOTYPES -DHAVE_OPENGL -DHKS_OPEN_GL \
###	    -DSUN_OGL_NO_VERTEX_MACROS -Winline -fpermissive \
###	    -fexceptions -DTYPENAME= -D_POSIX_SOURCE \
###	    -D_XOPEN_SOURCE -DSPECIALIZE -DNDEBUG -O3 -I%I") -Wuninitialized
###set (CMAKE_CXX_FLAGS_DEBUG "-DNDEBUG") ##-DNDEBUG -march=native
###set (CMAKE_CXX_FLAGS_RELEASE "-ffast-math") -O3  -O3

set (CMAKE_CXX_FLAGS "-ffast-math -DNDEBUG -libblas -fopenmp -fPIC -DPIC -w -O3 -I%I")
set (CMAKE_CXX_FLAGS_DEBUG "-DNDEBUG")
set (CMAKE_CXX_FLAGS_RELEASE "-ffast-math")

## target_include_directories(UserSubLib PRIVATE ${BLASROOT})
target_include_directories(UserSubLib PRIVATE ${EIGENROOT})
target_include_directories(UserSubLib PRIVATE ${DISROOT})
target_include_directories(UserSubLib PRIVATE ${ABQROOT} )

install(TARGETS UserSubLib DESTINATION ${DISROOT} )
