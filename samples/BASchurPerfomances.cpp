#include <iostream>
#include <vector>

#include "Problems.h"

#include <chrono>
#include <fstream>

#include "csNelson/DoubleSection.hpp"
#include "csNelson/DoubleSectionHessian.hpp"
#include "csNelson/EdgeBinary.hpp"
#include "csNelson/EdgeUnary.hpp"
#include "csNelson/GaussNewton.hpp"
#include "csNelson/LevenbergMarquardt.hpp"

#include "csCamera/Camera.hpp"
#include "csCamera/Utils.hpp"

#include "csLie/SE3_impl.hpp"
#include "csLie/Utils_impl.hpp"

constexpr int BR = 6;
constexpr int BC = 3;

template <int matUType, int matVType, int matWType>
class Edge;

template <int matUType, int matVType, int matWType>
class BA_Problem : public csNelson::DoubleSection<BA_Problem<matUType, matVType, matWType>, csLie::SE3d, Eigen::Vector3d, matUType,
                                                  matVType, matWType, double, 6, 3, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> {
  csCamera::Camerad _camera;
  std::vector<csLie::SE3d> _poses;
  std::vector<Eigen::Vector3d> _points;

  using DoubleSectionBase = csNelson::DoubleSection<BA_Problem<matUType, matVType, matWType>, csLie::SE3d, Eigen::Vector3d, matUType,
                                                    matVType, matWType, double, 6, 3, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>;

public:
  BA_Problem(int nc, int np, int ne)
      : _poses(nc)
      , _points(np)
      , _camera(100, 100, 50, 50, 100, 100) {
    const double radius_points = 1;
    for (int i = 0; i < _points.size(); i++) {
      _points[i] = Eigen::Vector3d::Random() * radius_points;
    }

    const double radius_cam = _camera.fx() * radius_points / (_camera.w() / 2.0);
    for (int i = 0; i < _poses.size(); i++) {
      Eigen::Vector3d cam_center = Eigen::Vector3d::Zero();
      while (cam_center.squaredNorm() < Eigen::NumTraits<double>::dummy_precision()) {
        cam_center.setRandom();
        cam_center.normalize();
      }

      Eigen::Vector3d dir_z = -cam_center;
      Eigen::Vector3d dir_r = Eigen::Vector3d::Zero();
      dir_r.setRandom();
      dir_r.normalize();
      while (1 - std::abs(dir_r.dot(dir_z)) < Eigen::NumTraits<double>::dummy_precision()) {
        dir_r.setRandom();
        dir_r.normalize();
      }
      Eigen::Vector3d dir_y = dir_z.cross(dir_r).normalized();
      Eigen::Vector3d dir_x = dir_y.cross(dir_z).normalized();

      _poses[i].linear().col(0) = dir_x;
      _poses[i].linear().col(1) = dir_y;
      _poses[i].linear().col(2) = dir_z;
      _poses[i].translation() = cam_center * radius_cam;

      _poses[i] = csLie::se3Exp(csLie::se3Log(_poses[i])).inverse();
    }

    this->reserveEdges(ne);
  }

  void addNoise() {

    for (int i = 1; i < _poses.size(); i++) {
      auto alg = csLie::SE3Algd();
      alg.head<3>() = Eigen::Vector3d() * 0.001;
      alg.tail<3>() = Eigen::Vector3d() * 0.001;
      _poses[i] = csLie::se3Exp(alg) * _poses[i];
    }

    for (int i = 1; i < _points.size(); i++) {
      _points[i] += (Eigen::Vector3d::Random() * 0.001);
    }
  }

  const csCamera::Camerad& camera() const {
    return _camera;
  }

  int numFixedParametersU() const override {
    return 1;
  }
  int numFixedParametersV() const override {
    return 1;
  }

  const csLie::SE3d& pose(int i) const {
    return _poses[i];
  }
  const Eigen::Vector3d& point(int i) const {
    return _points[i];
  }

  virtual const csLie::SE3d& parameterU(csNelson::NodeId i) const override {
    if (i.isVariable()) {
      return _poses[i.id() + 1];
    } else {
      assert(i.id() == 0);
      return _poses[i.id()];
    }
  }
  virtual csLie::SE3d& parameterU(csNelson::NodeId i) override {
    if (i.isVariable()) {
      return _poses[i.id() + 1];
    } else {
      assert(i.id() == 0);
      return _poses[i.id()];
    }
  }

  virtual const Eigen::Vector3d& parameterV(csNelson::NodeId i) const override {
    if (i.isVariable()) {
      return _points[i.id() + 1];
    } else {
      assert(i.id() == 0);
      return _points[i.id()];
    }
  }
  virtual Eigen::Vector3d& parameterV(csNelson::NodeId i) override {
    if (i.isVariable()) {
      return _points[i.id() + 1];
    } else {
      assert(i.id() == 0);
      return _points[i.id()];
    }
  }

  int numParametersV() const override {
    return _points.size() - 1;
  }

  int numParametersU() const override {
    return _poses.size() - 1;
  }

  void oplus(const typename DoubleSectionBase::HessianVectorsT& inc) {
    for (int i = 1; i < _points.size(); i++) {
      _points[i] += inc.bV().segment(i - 1);
    }
    for (int i = 1; i < _poses.size(); i++) {
      _poses[i] = csLie::se3Exp(csLie::SE3Algd(inc.bU().segment(this->user2internalIndexesU()(i - 1)))) * _poses[i];
    }
  }
};

template <int matUType, int matVType, int matWType>
class Edge : public BA_Problem<matUType, matVType, matWType>::template EdgeBinaryUV<Edge<matUType, matVType, matWType>> {
  Eigen::Vector2d _meas;

  Eigen::Vector3d _pC;
  Eigen::Vector2d _pZ1;
  Eigen::Vector2d _pImg;
  Eigen::Vector2d _err;

  Eigen::Matrix<double, 3, 6> _d_pC_d_gamma;
  Eigen::Matrix<double, 2, 3> _d_pZ1_d_pC;
  Eigen::Matrix<double, 2, 2> _d_pImg_d_pZ1;

  Eigen::Matrix<double, 2, 6> _d_err_d_gamma;
  Eigen::Matrix<double, 2, 3> _d_err_d_point;

public:
  Edge(const Eigen::Vector2d& meas)
      : _meas(meas) {}

  void update(bool hessians) override {
    const auto& pose = this->parameter_1();
    const auto& point = this->parameter_2();

    if (hessians) {
      this->_pC = pose * point;
      _d_pC_d_gamma = csLie::se3JacobianLeftZeroPerturbationPointTransformation(_pC);
      this->_pZ1 = csCamera::Utils::points3D_to_z1_jacob(_pC, _d_pZ1_d_pC);
      this->_pImg = this->section().camera().pointsZ1_to_image_jacobian(_pZ1, _d_pImg_d_pZ1);

      _d_err_d_gamma = _d_pImg_d_pZ1 * _d_pZ1_d_pC * _d_pC_d_gamma;
      _d_err_d_point = _d_pImg_d_pZ1 * _d_pZ1_d_pC * pose.linear();

      _err = this->_pImg - _meas;
      this->setChi2(_err.squaredNorm());
    } else {
      _err = this->section().camera().points3D_to_image(pose * point) - _meas;
      this->setChi2(_err.squaredNorm());
    }
  }

  template <class Derived1, class Derived2>
  void updateH11Block(Eigen::MatrixBase<Derived1>& H, Eigen::MatrixBase<Derived2>& v) {
    // wrt camera
    H.noalias() += _d_err_d_gamma.transpose() * _d_err_d_gamma;
    v += _d_err_d_gamma.transpose() * _err;
  }
  template <class Derived1, class Derived2>
  void updateH22Block(Eigen::MatrixBase<Derived1>& H, Eigen::MatrixBase<Derived2>& v) {
    H.noalias() += _d_err_d_point.transpose() * _d_err_d_point;
    v += _d_err_d_point.transpose() * _err;
  }
  template <class Derived1>
  void updateH12Block(Eigen::MatrixBase<Derived1>& H, bool transpose) {
    assert(transpose == false);
    H.noalias() += _d_err_d_gamma.transpose() * _d_err_d_point;
  }
};

template <int matUType, int matVType, int matWType>
void testFunction(const Problem& problem, bool AMDOrdering, int numThreads) {
  BA_Problem<matUType, matVType, matWType> ba_problem(problem.nCameras, problem.nPoints, problem.edges.size());
  ba_problem.parametersReady();
  ba_problem.reserveEdges(problem.edges.size());

  for (int i = 0; i < problem.edges.size(); i++) {
    const auto& e = problem.edges[i];
    const int e_i = e.first;
    const int e_j = e.second;
    csNelson::NodeId cam_i =
        csNelson::NodeId((e_i == 0 ? 0 : e_i - 1), (e_i == 0 ? csNelson::NodeType::Fixed : csNelson::NodeType::Variable));
    csNelson::NodeId point_i =
        csNelson::NodeId((e_j == 0 ? 0 : e_j - 1), (e_j == 0 ? csNelson::NodeType::Fixed : csNelson::NodeType::Variable));
    ba_problem.addEdge(cam_i, point_i,
                       new Edge<matUType, matVType, matWType>(
                           ba_problem.camera().points3D_to_image(ba_problem.pose(e.first) * ba_problem.point(e.second))));
  }

  if (AMDOrdering) {
    auto t0 = std::chrono::steady_clock::now();
    ba_problem.permuteAMD_SchurU();
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "PERMUTATION COMPUTE TIME " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  }

  ba_problem.structureReady();
  ba_problem.update(false);
  std::cout << "chi2 init " << ba_problem.hessian().chi2() << std::endl;

  ba_problem.addNoise();
  ba_problem.update(false);

  std::cout << "chi2 noise add " << ba_problem.hessian().chi2() << std::endl;

  using SolverAlgorithm = csNelson::GaussNewton<typename csNelson::SolverTraits<
      csNelson::solverCholeskySchurDiagBlockInverseWWtMult>::Solver<typename BA_Problem<matUType, matVType, matWType>::Hessian::Traits,
                                                                    csNelson::matrixWrapperSparse, csNelson::choleskyNaturalOrdering>>;
  SolverAlgorithm gn;

  ba_problem.settings().edgeEvalParallelSettings.setNumThreads(numThreads);
  ba_problem.settings().hessianUpdateParallelSettings.setNumThreads(numThreads);
  gn.solverSettings().Vinv.blockInversion.setNumThreads(numThreads);
  gn.solverSettings().Vinv.rightVectorMult.setNumThreads(numThreads);
  gn.solverSettings().WtX.setNumThreads(numThreads);
  gn.solverSettings().WVinv.multiplication.setNumThreads(numThreads);
  gn.solverSettings().WVinv.rightVectorMult.setNumThreads(numThreads);
  gn.solverSettings().WVinvWt.setNumThreads(numThreads);
  gn.settings().minNumIt = 5;

  auto t0 = std::chrono::steady_clock::now();
  auto tc = gn.solve(ba_problem);
  auto t1 = std::chrono::steady_clock::now();

  std::cout << "ELAPSED " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
  std::cout << "stats " << gn.stats().toString() << std::endl;
  std::cout << gn.timingStats().toString() << std::endl;
  std::cout << gn.solver().timingStats().toString() << std::endl;
}

void test(const Problem& problem, bool amdOrdering, int numThreads) {
  std::cout << "----------------------------------------------------" << std::endl
            << "U BlockDiagonal, V BlockDiagonal, W BlockSparse" << std::endl;
  testFunction<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>(problem, amdOrdering,
                                                                                                       numThreads);

  std::cout << "----------------------------------------------------" << std::endl
            << "U SparseCoeffBlockDiagonal, V SparseCoeffBlockDiagonal, W BlockCoeffSparse" << std::endl;
  testFunction<csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::BlockCoeffSparse>(
      problem, amdOrdering, numThreads);
  std::cout << "----------------------------------------------------" << std::endl << std::endl;

  std::cout << "----------------------------------------------------" << std::endl
            << "U BlockDiagonal, V BlockDiagonal, W BlockCoeffSparse" << std::endl;
  testFunction<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>(problem, amdOrdering,
                                                                                                            numThreads);
  std::cout << "----------------------------------------------------" << std::endl << std::endl;

  std::cout << "----------------------------------------------------" << std::endl
            << "U SparseCoeffBlockDiagonal, V SparseCoeffBlockDiagonal, W BlockSparse" << std::endl;
  testFunction<csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::BlockSparse>(
      problem, amdOrdering, numThreads);
  std::cout << "----------------------------------------------------" << std::endl << std::endl;
}

int main(int argc, char* argv[]) {

  const Problem* problemPtr = &ProblemCollections::ladybug_49_7776;
  std::unique_ptr<Problem> newProblem;
  if (argc == 2) {
    newProblem.reset(new Problem());
    bool ok = newProblem->load(argv[1]);
    if (!ok) {
      std::cerr << "error reading " << argv[1] << std::endl;
      std::exit(-1);
    } else {
      problemPtr = newProblem.get();
    }
  }
  const Problem& problem = *problemPtr;

  // int numThreads = csNelson::ParallelExecSettings::maxSupportedThreads() / 2;
  int numThreads = 1;

  std::cout << "********* AMD ********* " << std::endl;
  test(problem, true, numThreads);
  std::cout << "*********************** " << std::endl << std::endl << std::endl;
  std::cout << "********* NAT ********* " << std::endl;
  test(problem, false, numThreads);
  std::cout << "*********************** " << std::endl << std::endl << std::endl;

  return 0;
}
