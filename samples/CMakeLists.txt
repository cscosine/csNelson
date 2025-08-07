project(csNelsonSamples LANGUAGES CXX)

#---------------------------------------------------------------
# prepare SamplesProblemLib
cs_add_library(csNelsonSamplesProblem STATIC)
cs_add_library(csNelsonSamplesProblem::csNelsonSamplesProblem ALIAS csNelsonSamplesProblem)

cs_target_sources(csNelsonSamplesProblem
	PRIVATE_INSTALL
    Problems.hpp
	PRIVATE
    Problems.cpp
)
cs_target_link_libraries(csNelsonSamplesProblem PUBLIC csNelsonSol::csNelsonSol_csCMakeBaseInclude)

cs_install_target(csNelsonSamplesProblem csNelsonSamplesProblem DEFAULT)

cs_finalize_target(csNelsonSamplesProblem)

#---------------------------------------------------------------
# prepare a set of samples dependent on the above library
set(SAMPLES_LIST
  WWtMultPerformances.cpp
  WtXMultPerformances.cpp
  WVinvMultPerformances.cpp
  VinvMultPerformances.cpp
  TestSetUnorderedSetPerformances.cpp
)

foreach(s ${SAMPLES_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	cs_add_executable(csNelsonSamples${sn})
	cs_target_sources(csNelsonSamples${sn} PRIVATE ${s})
	cs_target_link_libraries(csNelsonSamples${sn} PRIVATE csNelson::csNelson csNelsonSamplesProblem::csNelsonSamplesProblem)
	cs_install_target(csNelsonSamples${sn} csNelsonSamples${sn})
	cs_finalize_target(csNelsonSamples${sn})
endforeach()


#---------------------------------------------------------------
# prepare a second set of samples, with some additional dependencies from real problems (lie and camera)
set(SAMPLES_LIST
  BASchurPerfomances.cpp
)

foreach(s ${SAMPLES_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	cs_add_executable(csNelsonSamples${sn})
	cs_find_package(csNelsonSamples${sn} PUBLIC csLie COMPONENTS csLie REQUIRED)
	cs_find_package(csNelsonSamples${sn} PUBLIC csCamera COMPONENTS csCamera REQUIRED)
	cs_target_sources(csNelsonSamples${sn} PRIVATE ${s})
	cs_target_link_libraries(csNelsonSamples${sn} PRIVATE csNelson::csNelson csNelsonSamplesProblem::csNelsonSamplesProblem csLie::csLie csCamera::csCamera)
	cs_install_target(csNelsonSamples${sn} csNelsonSamples${sn})
	cs_finalize_target(csNelsonSamples${sn})
endforeach()

#---------------------------------------------------------------
# prepare a set of independent samples
set(SAMPLES_LIST
  nelsonPlayground.cpp
)

foreach(s ${SAMPLES_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	cs_add_executable(csNelsonSamples${sn})
	cs_target_sources(csNelsonSamples${sn} PRIVATE ${s})
	cs_target_link_libraries(csNelsonSamples${sn} PRIVATE csNelson::csNelson)
	cs_install_target(csNelsonSamples${sn} csNelsonSamples${sn})
	cs_finalize_target(csNelsonSamples${sn})
endforeach()
