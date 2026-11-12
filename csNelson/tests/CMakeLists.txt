project(csNelsonTest  LANGUAGES CXX)

include(CTest)

set(TESTS_LIST
  testSingleSectionHessian.cpp
  testDoubleSectionHessian.cpp
  testSingleSection.cpp
  testSingleSectionPermutation.cpp
  testMatrixSparseDenseWrapper.cpp
  testMatrixWWtMultiplier.cpp
  testMatrixWVinvMultiplier.cpp
  testMatrixWtXMultiplier.cpp
  testMatrixVinv.cpp
  testDoubleSectionSchurDiagBlockWWtFFFF.cpp
  testEdgeNary.cpp
)

foreach(s ${TESTS_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	cs_add_executable(csNelsonTest${sn})
  cs_find_package(csNelsonTest${sn} PUBLIC Catch2 3 REQUIRED CONFIG)
  include(Catch)
	cs_target_sources(csNelsonTest${sn} PRIVATE ${s})
	cs_target_link_libraries(csNelsonTest${sn} PRIVATE csNelson::csNelson Catch2::Catch2WithMain )
  catch_discover_tests(csNelsonTest${sn} DISCOVERY_MODE "PRE_TEST")
  cs_install_target(csNelsonTest${sn} COMPONENT csNelsonTest${sn})
endforeach()

######################################################################
OPTION(csNelson_COMPILE_LONG_TESTS "Compile All Tests, including very long compile time ones" ON)

if(csNelson_COMPILE_LONG_TESTS)
  set(DOUBLESECTIONTESTCOMMON testDoubleSectionCommon.h testDoubleSectionCommon.hpp)

  set(TESTS_LIST
    testDoublesectionFFFF.cpp
    testDoublesectionFFFD.cpp
    testDoublesectionFFDF.cpp
    testDoublesectionFFDD.cpp
    #------------------------
    testDoublesectionFDFF.cpp
    testDoublesectionFDFD.cpp
    testDoublesectionFDDF.cpp
    testDoublesectionFDDD.cpp
    #------------------------
    testDoublesectionFVFF.cpp
    testDoublesectionFVFD.cpp
    testDoublesectionFVDF.cpp
    testDoublesectionFVDD.cpp
    #------------------------
    #------------------------
    testDoublesectionDFFF.cpp
    testDoublesectionDFFD.cpp
    testDoublesectionDFDF.cpp
    testDoublesectionDFDD.cpp
    #------------------------
    testDoublesectionDDFF.cpp
    testDoublesectionDDFD.cpp
    testDoublesectionDDDF.cpp
    testDoublesectionDDDD.cpp
    #------------------------
    testDoublesectionDVFF.cpp
    testDoublesectionDVFD.cpp
    testDoublesectionDVDF.cpp
    testDoublesectionDVDD.cpp
    #------------------------
    #------------------------
    testDoublesectionVFFF.cpp
    testDoublesectionVFFD.cpp
    testDoublesectionVFDF.cpp
    testDoublesectionVFDD.cpp
    #------------------------
    testDoublesectionVDFF.cpp
    testDoublesectionVDFD.cpp
    testDoublesectionVDDF.cpp
    testDoublesectionVDDD.cpp
    #------------------------
    testDoublesectionVVFF.cpp
    testDoublesectionVVFD.cpp
    testDoublesectionVVDF.cpp
    testDoublesectionVVDD.cpp
  )

  foreach(s ${TESTS_LIST})
    get_filename_component(sn "${s}" NAME_WE)
    cs_add_executable(csNelsonTest${sn})
    cs_find_package(csNelsonTest${sn} PUBLIC Catch2 3 REQUIRED CONFIG)
    include(Catch)
    cs_target_sources(csNelsonTest${sn} PRIVATE ${s} ${DOUBLESECTIONTESTCOMMON})
    cs_target_link_libraries(csNelsonTest${sn} PRIVATE csNelson::csNelson Catch2::Catch2WithMain )
    catch_discover_tests(csNelsonTest${sn} DISCOVERY_MODE "PRE_TEST")
    cs_install_target(csNelsonTest${sn} COMPONENT csNelsonTest${sn})
  endforeach()


#--------------------------------------------------------------------------------------------


  set(DOUBLESECTIONSCHURTESTCOMMON testDoubleSectionSchurCommon.h testDoubleSectionSchurCommon.hpp)


  set(TESTS_LIST
    testDoubleSectionSchurDiagBlockFFFF.cpp
    testDoubleSectionSchurDiagBlockFFFD.cpp
    testDoubleSectionSchurDiagBlockFFDF.cpp
    testDoubleSectionSchurDiagBlockFFDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockFDFF.cpp
    testDoubleSectionSchurDiagBlockFDFD.cpp
    testDoubleSectionSchurDiagBlockFDDF.cpp
    testDoubleSectionSchurDiagBlockFDDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockFVFF.cpp
    testDoubleSectionSchurDiagBlockFVFD.cpp
    testDoubleSectionSchurDiagBlockFVDF.cpp
    testDoubleSectionSchurDiagBlockFVDD.cpp
    #--------------------------------------
    #--------------------------------------
    testDoubleSectionSchurDiagBlockDFFF.cpp
    testDoubleSectionSchurDiagBlockDFFD.cpp
    testDoubleSectionSchurDiagBlockDFDF.cpp
    testDoubleSectionSchurDiagBlockDFDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockDDFF.cpp
    testDoubleSectionSchurDiagBlockDDFD.cpp
    testDoubleSectionSchurDiagBlockDDDF.cpp
    testDoubleSectionSchurDiagBlockDDDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockDVFF.cpp
    testDoubleSectionSchurDiagBlockDVFD.cpp
    testDoubleSectionSchurDiagBlockDVDF.cpp
    testDoubleSectionSchurDiagBlockDVDD.cpp
    #--------------------------------------
    #--------------------------------------
    testDoubleSectionSchurDiagBlockVFFF.cpp
    testDoubleSectionSchurDiagBlockVFFD.cpp
    testDoubleSectionSchurDiagBlockVFDF.cpp
    testDoubleSectionSchurDiagBlockVFDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockVDFF.cpp
    testDoubleSectionSchurDiagBlockVDFD.cpp
    testDoubleSectionSchurDiagBlockVDDF.cpp
    testDoubleSectionSchurDiagBlockVDDD.cpp
    #--------------------------------------
    testDoubleSectionSchurDiagBlockVVFF.cpp
    testDoubleSectionSchurDiagBlockVVFD.cpp
    testDoubleSectionSchurDiagBlockVVDF.cpp
    testDoubleSectionSchurDiagBlockVVDD.cpp
  )

  foreach(s ${TESTS_LIST})
    get_filename_component(sn "${s}" NAME_WE)
    cs_add_executable(csNelsonTest${sn})
    cs_find_package(csNelsonTest${sn} PUBLIC Catch2 3 REQUIRED CONFIG)
    include(Catch)
    cs_target_sources(csNelsonTest${sn} PRIVATE ${s} ${DOUBLESECTIONSCHURTESTCOMMON})
    cs_target_link_libraries(csNelsonTest${sn} PRIVATE csNelson::csNelson Catch2::Catch2WithMain)
    catch_discover_tests(csNelsonTest${sn} DISCOVERY_MODE "PRE_TEST")
    cs_install_target(csNelsonTest${sn} COMPONENT csNelsonTest${sn})
  endforeach()
endif()

######################################################################

set(TESTS_LIST
  testGaussNewton2DScansAlign.cpp
  testGaussNewton2DScansAlignNary.cpp
)

foreach(s ${TESTS_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	cs_add_executable(csNelsonTest${sn})
  cs_find_package(csNelsonTest${sn} PUBLIC Catch2 3 REQUIRED CONFIG)
  include(Catch)
  cs_find_package(csNelsonTest${sn} PUBLIC csLie COMPONENTS csLie REQUIRED CONFIG)
	cs_target_sources(csNelsonTest${sn} PRIVATE ${s})
	cs_target_link_libraries(csNelsonTest${sn} PRIVATE csNelson::csNelson Catch2::Catch2WithMain csLie::csLie)
  catch_discover_tests(csNelsonTest${sn} DISCOVERY_MODE "PRE_TEST")
  cs_install_target(csNelsonTest${sn} COMPONENT csNelsonTest${sn})
endforeach()
