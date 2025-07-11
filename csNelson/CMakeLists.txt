project(csNelson LANGUAGES CXX)


# prepare custom project asserts, note: use standard find_package
find_package(csCore REQUIRED COMPONENTS assert CONFIG)

list(APPEND CMAKE_MODULE_PATH "${csCore_BASE_DIR}/cmake/csCore/assert/")
include(csAssert)

# set defaults for Paranoid additional configuration
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT RELEASE Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT RELWITHDEBINFO Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT DEBUG Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT PARANOID Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT SAFE_MATH Paranoid ON)

csAssertGenerateProject(csNelson csNelson) #generate target csNelson::csNelsonAssert
list(POP_BACK CMAKE_MODULE_PATH)

# prepare standard library
cs_add_library(csNelson STATIC)
cs_add_library(csNelson::csNelson ALIAS csNelson)

cs_find_package(csNelson PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csNelson PUBLIC csCore COMPONENTS core eigenUtils REQUIRED CONFIG)

message(STATUS "Looking for OpenMP. If this fail, for clang try `sudo apt install libomp-dev`")
cs_find_package(csNelson PUBLIC OpenMP REQUIRED)

cs_find_package(csNelson PUBLIC csBlockMatrix COMPONENTS csBlockMatrix REQUIRED)

cs_target_sources(csNelson
	PRIVATE_INSTALL
    Global.h
    BaseSection.h
    SingleSectionHessian.h
    SingleSectionHessian.hpp
    DoubleSectionHessianMatrices.h
    DoubleSectionHessianMatrices.hpp
    DoubleSectionHessian.h
    DoubleSectionHessian.hpp
    SingleSection.h
    SingleSection.hpp
    DoubleSection.h
    DoubleSection.hpp
    EdgeInterface.h
    EdgeSectionBase.h
    EdgeSectionBase.hpp
    EdgeNary.h
    EdgeNary.hpp
    EdgeUnary.h
    EdgeUnary.hpp
    EdgeBinary.h
    EdgeBinary.hpp
    MatrixWrapperTraits.h
    MatrixDenseWrapper.h
    MatrixDenseWrapper.hpp
    MatrixSparseWrapper.h
    MatrixSparseWrapper.hpp
    MatrixDiagInv.h
    MatrixDiagInv.hpp
    MatrixWWtMultiplier.h
    MatrixWWtMultiplier.hpp
    MatrixWVinvMultiplier.h
    MatrixWVinvMultiplier.hpp
    MatrixWtXMultiplier.h
    MatrixWtXMultiplier.hpp
    GaussNewton.h
    GaussNewton.hpp
    SolverCholeskyDense.h
    SolverCholeskyDense.hpp
    SolverCholeskySchur.h
    SolverCholeskySchur.hpp
    SolverCholeskySparse.h
    SolverCholeskySparse.hpp
    SolverDiagonalBlocksInverseSchur.h
    SolverDiagonalBlocksInverseSchur.hpp
    SolverDiagonalBlocksInverseWWtMultSchur.h
    SolverDiagonalBlocksInverseWWtMultSchur.hpp
    SolverTraitsBase.h
    SolverTraits.h
    ParallelExecHelper.h
    LevenbergMarquardt.h
    LevenbergMarquardt.hpp
    OrderingTraits.h
	PRIVATE
    Global.cpp
    BaseSection.cpp
    SingleSectionHessian.cpp
    DoubleSectionHessian.cpp
    SingleSection.cpp
    DoubleSection.cpp
    EdgeSectionBase.cpp
    EdgeInterface.cpp
    EdgeUnary.cpp
    EdgeNary.cpp
    EdgeBinary.cpp
    MatrixDenseWrapper.cpp
    MatrixSparseWrapper.cpp
    MatrixDiagInv.cpp
    MatrixWWtMultiplier.cpp
    MatrixWVinvMultiplier.cpp
    MatrixWtXMultiplier.cpp
    GaussNewton.cpp
    SolverCholeskyDense.cpp
    SolverCholeskySchur.cpp
    SolverCholeskySparse.cpp
    SolverDiagonalBlocksInverseSchur.cpp
    SolverDiagonalBlocksInverseWWtMultSchur.cpp
    ParallelExecHelper.cpp
    LevenbergMarquardt.cpp
)

cs_target_link_libraries(csNelson
	PUBLIC
		Eigen3::Eigen
    	csCore::core
      csBlockMatrix::csBlockMatrix
		csCore::eigenUtils
		csNelson::csNelsonAssert
		csNelsonLib::csNelsonLib_csCMakeBaseInclude
)

# TODO(sc) re-enable tests
# add_subdirectory(tests)

cs_install_target(csNelson csNelson DEFAULT)

cs_finalize_target(csNelson)
