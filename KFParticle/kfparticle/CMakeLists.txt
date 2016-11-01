cmake_minimum_required(VERSION 2.6.2)
project(CA)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Vc REQUIRED)
include_directories(${VC_INCLUDE_DIR})

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebug RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(Root REQUIRED)

if(NUM_THREADS)
   add_definitions(-DNUM_THREADS=${NUM_THREADS})
endif(NUM_THREADS)
if(VC_NO_GATHER_TRICKS)
   add_definitions(-DVC_NO_GATHER_TRICKS)
endif(VC_NO_GATHER_TRICKS)

macro(add_target_property _target _prop _value)
   get_target_property(_oldprop "${_target}" ${_prop})
   if(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
   else(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
   endif(NOT _oldprop)
endmacro(add_target_property)

add_definitions(-DDO_TPCCATRACKER_EFF_PERFORMANCE) # when comment it disable ROOT in 'cmake -i' as well
add_definitions(-DHomogeneousField)
add_definitions(-DUSE_TIMERS)
add_definitions(-DKFParticleStandalone)
    
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse -msse2")

if(CMAKE_COMPILER_IS_GNUCXX)
   ##################################################################################################
   #                                              GCC                                               #
   ##################################################################################################
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -W -Wall -Wswitch -Wformat -Wchar-subscripts -Wno-parentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wno-variadic-macros -Wshadow -fno-threadsafe-statics -fopenmp")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -W -Wall -Wswitch -Wformat -Wchar-subscripts -Wno-parentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wno-variadic-macros -Wshadow -fno-threadsafe-statics -fopenmp")

   if(ENABLE_TBB)
      set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -DTBB_DO_ASSERT=1")
   endif(ENABLE_TBB)
   set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}")
   set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")
   if(ENABLE_TBB)
      set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -DTBB_DO_ASSERT=1")
   endif(ENABLE_TBB)
   set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}")
   set(CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -O3")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -O3")

   
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")

   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      add_definitions(-DNO_OPTIMIZATION)
      add_definitions(-g)
   elseif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
      add_definitions(-DNO_OPTIMIZATION)
   endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
   if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
      set(ENABLE_STRICT_ALIASING true CACHE BOOL "Enables strict aliasing rules for more aggressive optimizations")
      if(NOT ENABLE_STRICT_ALIASING)
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")
         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing")
      endif(NOT ENABLE_STRICT_ALIASING)
      set(DISABLE_DEBUG true CACHE BOOL "Disables debugging code, like assertions")
      if(DISABLE_DEBUG)
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -DNODEBUG")
         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG -DNODEBUG")
      endif(DISABLE_DEBUG)
   endif(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")

   set(USE_SSSE3 ${SSSE3_FOUND} CACHE BOOL "Use SSE3 and SSSE3. If (S)SSE3 instructions are not enabled they will be emulated.")
   set(USE_SSE4_1 ${SSE4_1_FOUND} CACHE BOOL "Use SSE4.1. If SSE4.1 instructions are not enabled they will be emulated.")
   if(USE_SSSE3)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -mssse3")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -mssse3")
   endif(USE_SSSE3)
   if(USE_SSE4_1)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
   endif(USE_SSE4_1)
#   set(USE_AVX ${AVX_FOUND} CACHE BOOL "Use AVX. If AVX instructions are not enabled they will be emulated.")
   if(AVX_FOUND)
      message(STATUS "AVX is found (disabled by default).")
      set(USE_AVX FALSE CACHE BOOL "Use AVX. (AVX is supported and can be anabled).")
   elseif(AVX_FOUND)
      set(USE_AVX FALSE CACHE BOOL "Use AVX. (AVX is not supported).")
   endif(AVX_FOUND)
   if(USE_AVX)
      message(STATUS "AVX will be used.")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
   endif(USE_AVX)
elseif(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")
   ##################################################################################################
   #                                          Intel Compiler                                        #
   ##################################################################################################
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -W -Wall -Wformat -Wno-parentheses -Wmultichar -Wtrigraphs -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wno-variadic-macros -Wshadow -std=c++11 -openmp")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -W -Wall -Wformat -Wno-parentheses -Wmultichar -Wtrigraphs -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wno-variadic-macros -Wshadow -std=c++11 -openmp")

   # add_definitions(-Wno-pragmas) #   root has unknown for ICC pragmas
  
   if(ENABLE_TBB)   
      set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -DTBB_DO_ASSERT=1")
   endif(ENABLE_TBB)

   set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -DNODEBUG -O3")
   if(ENABLE_TBB)
      set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -DTBB_DO_ASSERT=1")
   endif(ENABLE_TBB)
   set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE} -O3")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -DNODEBUG -O3")

   set(ALIAS_FLAGS "-no-ansi-alias")
   if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
      set(ENABLE_STRICT_ALIASING true CACHE BOOL "Enables strict aliasing rules for more aggressive optimizations")
      if(ENABLE_STRICT_ALIASING)
         set(ALIAS_FLAGS "-ansi-alias")
      endif(ENABLE_STRICT_ALIASING)
   endif(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")

   # per default icc is not IEEE compliant, but we need that for verification
   set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${ALIAS_FLAGS} -w1 -fp-model precise")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ALIAS_FLAGS} -w1 -fp-model precise")

   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      add_definitions(-DNO_OPTIMIZATION)
      add_definitions(-g)
   elseif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
      add_definitions(-DNO_OPTIMIZATION)
   endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

   # set(USE_SSSE3 ${SSSE3_FOUND} CACHE BOOL "Use SSE3 and SSSE3. If (S)SSE3 instructions are not enabled they will be emulated.")
   # set(USE_SSE4_1 ${SSE4_1_FOUND} CACHE BOOL "Use SSE4.1. If SSE4.1 instructions are not enabled they will be emulated.")
   # if(USE_SSSE3)
   #    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -mssse3")
   #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -mssse3")
   # endif(USE_SSSE3)
   # if(USE_SSE4_1)
   #    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
   #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
   # endif(USE_SSE4_1)
   if(COMPILE_FOR_MIC)
   else(COMPILE_FOR_MIC)
    if(AVX_FOUND)
        message(STATUS "AVX is found (enabled by default).")
        set(USE_AVX TRUE CACHE BOOL "Use AVX. (AVX is supported and can be anabled).")
    elseif(AVX_FOUND)
        set(USE_AVX FALSE CACHE BOOL "Use AVX. (AVX is not supported).")
    endif(AVX_FOUND)
    if(USE_AVX)
        message(STATUS "AVX will be used.")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xAVX")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xAVX")
    endif(USE_AVX)
  endif(COMPILE_FOR_MIC)
endif(CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
   message(STATUS "WARNING! It seems you are compiling without optimization. Please set CMAKE_BUILD_TYPE.")
endif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")

set(FAST_MATH FALSE CACHE BOOL "Use floating point math shortcuts. Might not be IEEE compliant!")
if(FAST_MATH)
   if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
   elseif(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model fast=2")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model fast=2")
   elseif(MSVC)
   endif(CMAKE_COMPILER_IS_GNUCXX)
endif(FAST_MATH)

set(PROFILING FALSE CACHE BOOL "Compile the code for profiling with gprof. Can slow down the code.")
if(PROFILING)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(PROFILING)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/KFParticle
  ${CMAKE_CURRENT_SOURCE_DIR}/KFParticlePerformance
)

if(ROOT_FOUND)
   include_directories(${ROOT_INCLUDE_DIR})
endif(ROOT_FOUND)

set(KFParticleCode
   KFParticle/KFParticleTopoReconstructor.cxx
   KFParticle/KFVertex.cxx	
   KFParticle/KFPTrack.cxx
   KFParticle/KFPTrackVector.cxx
   KFParticle/KFPVertex.cxx
   KFParticle/KFParticlePVReconstructor.cxx
   KFParticle/KFParticleDatabase.cxx
   KFParticle/KFParticleBase.cxx
   KFParticle/KFParticleBaseSIMD.cxx
   KFParticle/KFParticle.cxx
   KFParticle/KFParticleSIMD.cxx
   KFParticle/KFParticleFinder.cxx
)
set(KFPPerformanceCode
   KFParticlePerformance/KFMCVertex.cxx
   KFParticlePerformance/KFParticlePerformanceBase.cxx
   KFParticlePerformance/KFTopoPerformance.cxx
   KFParticlePerformance/KFMCParticle.cxx
)

set(USE_SIMD TRUE CACHE BOOL "When compiling on x86 use SSE Vector class instead of Simple implementation")
if(USE_SIMD)
# if(NOT COMPILE_FOR_MIC)
#    add_definitions(-DVC_IMPL=SSE)
# endif(NOT COMPILE_FOR_MIC)
else(USE_SIMD)
   add_definitions(-DVC_IMPL=Scalar)
endif(USE_SIMD)
set(CATrackerCode ${CATrackerCode} code/CATracker/AliHLTTPCCASliceDataVector.cxx)

#########################################################################################

#add_executable(convertToBinary convertToBinary.cpp code/AliHLTTPCCAMCTrack.cxx)
#target_link_libraries(convertToBinary ${ROOT_LIBS} ${VC_LIBRARIES})

#add_executable(benchmark benchmark.cpp)


add_library(KFParticle SHARED ${KFParticleCode} ${KFPPerformanceCode})
target_link_libraries(KFParticle ${VC_LIBRARIES} ${ROOT_LIBS})
if(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")
  add_target_property(KFParticle COMPILE_FLAGS "-inline-forceinline")
endif(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")

# install dir is not interesting for now
mark_as_advanced(CMAKE_INSTALL_PREFIX)
