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
    Global.hpp
    BaseSection.hpp
    SingleSectionHessian.hpp
    SingleSectionHessian_impl.hpp
    DoubleSectionHessianMatrices.hpp
    DoubleSectionHessianMatrices_impl.hpp
    DoubleSectionHessian.hpp
    DoubleSectionHessian_impl.hpp
    SingleSection.hpp
    SingleSection_impl.hpp
    DoubleSection.hpp
    DoubleSection_impl.hpp
    EdgeInterface.hpp
    EdgeSectionBase.hpp
    EdgeSectionBase_impl.hpp
    EdgeNary.hpp
    EdgeNary_impl.hpp
    EdgeUnary.hpp
    EdgeUnary_impl.hpp
    EdgeBinary.hpp
    EdgeBinary_impl.hpp
    MatrixWrapperTraits.hpp
    MatrixDenseWrapper.hpp
    MatrixDenseWrapper_impl.hpp
    MatrixSparseWrapper.hpp
    MatrixSparseWrapper_impl.hpp
    MatrixDiagInv.hpp
    MatrixDiagInv_impl.hpp
    MatrixWWtMultiplier.hpp
    MatrixWWtMultiplier_impl.hpp
    MatrixWVinvMultiplier.hpp
    MatrixWVinvMultiplier_impl.hpp
    MatrixWtXMultiplier.hpp
    MatrixWtXMultiplier_impl.hpp
    GaussNewton.hpp
    GaussNewton_impl.hpp
    SolverCholeskyDense.hpp
    SolverCholeskyDense_impl.hpp
    SolverCholeskySchur.hpp
    SolverCholeskySchur_impl.hpp
    SolverCholeskySparse.hpp
    SolverCholeskySparse_impl.hpp
    SolverDiagonalBlocksInverseSchur.hpp
    SolverDiagonalBlocksInverseSchur_impl.hpp
    SolverDiagonalBlocksInverseWWtMultSchur.hpp
    SolverDiagonalBlocksInverseWWtMultSchur_impl.hpp
    SolverTraitsBase.hpp
    SolverTraits.hpp
    ParallelExecHelper.hpp
    LevenbergMarquardt.hpp
    LevenbergMarquardt_impl.hpp
    OrderingTraits.hpp
	PRIVATE
    BaseSection.cpp
    EdgeInterface.cpp
    EdgeUnary.cpp
    EdgeBinary.cpp
    GaussNewton.cpp
    SolverCholeskySchur.cpp
    SolverDiagonalBlocksInverseSchur.cpp
    SolverDiagonalBlocksInverseWWtMultSchur.cpp
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
