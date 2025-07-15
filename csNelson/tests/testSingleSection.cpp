#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csBlockMatrix/DenseMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparsityPattern.hpp>

#include "csNelson/EdgeBinary.hpp"
#include "csNelson/EdgeUnary.hpp"
#include "csNelson/SingleSection.hpp"
#include "csNelson/SingleSectionHessian.hpp"

#include <array>
#include <iostream>

struct Point {
  Eigen::Vector3d point;
};

static constexpr int secSizeFix = 3;
static constexpr int numBlocks = 10;

template <class Section>
class EdgeUnaryTest : public csNelson::EdgeUnarySectionBaseCRPT<Section, typename Section::EdgeUnaryAdapter, EdgeUnaryTest<Section>> {
  int _parId;

public:
  EdgeUnaryTest(int parId)
      : _parId(parId) {}

  void update(bool hessians) override {
    if (this->parId().isVariable()) {
      REQUIRE(this->parId().id() == _parId);
      REQUIRE(this->HUid() >= 0);
    } else {
      REQUIRE(this->parId().id() == _parId);
    }

    const auto& par = this->parameter();
  }

  template <class Derived1, class Derived2>
  void updateHBlock(Eigen::MatrixBase<Derived1>& H, Eigen::MatrixBase<Derived2>& v) {
    REQUIRE(this->parId().isVariable());
    std::cout << "EdgeUnaryTest::updateHBlock " << this->parId().id() << "," << this->parId().id() << std::endl;
    H.setConstant(1);
    v.setConstant(1);
  }
};

template <class Section>
class EdgeBinaryTest
    : public csNelson::EdgeBinarySectionBaseCRPT<Section, typename Section::EdgeBinaryAdapter, EdgeBinaryTest<Section>> {
  int _par1Id, _par2Id;

public:
  EdgeBinaryTest(int par1Id, int par2Id)
      : _par1Id(par1Id)
      , _par2Id(par2Id) {}

  void update(bool hessians) override {
    REQUIRE(this->par_1_Id().id() == _par1Id);
    REQUIRE(this->par_2_Id().id() == _par2Id);
    if (this->par_1_Id().isVariable()) {
      REQUIRE(this->H_11_Uid() >= 0);
    }
    if (this->par_1_Id().isVariable() && this->par_2_Id().isVariable()) {
      REQUIRE(this->H_12_Uid() >= 0);
    }
    if (this->par_2_Id().isVariable()) {
      REQUIRE(this->H_22_Uid() >= 0);
    }

    const auto& p1 = this->parameter_1();
    const auto& p2 = this->parameter_2();
  }

  template <class Derived1, class Derived2>
  void updateH11Block(Eigen::MatrixBase<Derived1>& H, Eigen::MatrixBase<Derived2>& b) {
    REQUIRE(this->par_1_Id().isVariable());
    std::cout << "EdgeBinaryTest::updateHBlock[11] " << this->par_1_Id().id() << "," << this->par_1_Id().id() << std::endl;
    H.setConstant(11);
    b.setConstant(11);
  }
  template <class Derived>
  void updateH12Block(Eigen::MatrixBase<Derived>& b, bool transpose) {
    assert(transpose == false);
    REQUIRE(this->par_1_Id().isVariable());
    REQUIRE(this->par_2_Id().isVariable());
    std::cout << "EdgeBinaryTest::updateHBlock[12] " << this->par_1_Id().id() << "," << this->par_2_Id().id() << std::endl;
    b.setConstant(12);
  }
  template <class Derived1, class Derived2>
  void updateH22Block(Eigen::MatrixBase<Derived1>& H, Eigen::MatrixBase<Derived2>& b) {
    REQUIRE(this->par_2_Id().isVariable());
    std::cout << "EdgeBinaryTest::updateHBlock[22] " << this->par_2_Id().id() << "," << this->par_2_Id().id() << std::endl;
    H.setConstant(22);
    b.setConstant(22);
  }
};

template <int matType>
class PointsSectionFF : public csNelson::SingleSection<PointsSectionFF<matType>, Point, matType, double, secSizeFix, numBlocks> {
  std::array<Point, numBlocks> _points;
  Point _fixedPoint;

public:
  PointsSectionFF() {
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionFF_BlockDense = PointsSectionFF<csBlockMatrix::BlockDense>;
using PointsSectionFF_BlockDiagonal = PointsSectionFF<csBlockMatrix::BlockDiagonal>;
using PointsSectionFF_BlockSparse = PointsSectionFF<csBlockMatrix::BlockSparse>;
using PointsSectionFF_BlockCoeffSparse = PointsSectionFF<csBlockMatrix::BlockCoeffSparse>;

//-----------------------------------------------------------------------------------------------------------------

template <int matType>
class PointsSectionFD
    : public csNelson::SingleSection<PointsSectionFD<matType>, Point, matType, double, secSizeFix, csBlockMatrix::Dynamic> {
  std::vector<Point> _points;
  Point _fixedPoint;

public:
  PointsSectionFD() {
    _points.resize(numBlocks);
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  int numParameters() const override {
    return _points.size();
  }

  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionFD_BlockDense = PointsSectionFD<csBlockMatrix::BlockDense>;
using PointsSectionFD_BlockDiagonal = PointsSectionFD<csBlockMatrix::BlockDiagonal>;
using PointsSectionFD_BlockSparse = PointsSectionFD<csBlockMatrix::BlockSparse>;
using PointsSectionFD_BlockCoeffSparse = PointsSectionFD<csBlockMatrix::BlockCoeffSparse>;

//-----------------------------------------------------------------------------------------------------------------

template <int matType>
class PointsSectionDF
    : public csNelson::SingleSection<PointsSectionDF<matType>, Point, matType, double, csBlockMatrix::Dynamic, numBlocks> {
  std::array<Point, numBlocks> _points;
  Point _fixedPoint;

public:
  PointsSectionDF() {
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  int parameterSize() const override {
    return secSizeFix;
  }
  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionDF_BlockDense = PointsSectionDF<csBlockMatrix::BlockDense>;
using PointsSectionDF_BlockDiagonal = PointsSectionDF<csBlockMatrix::BlockDiagonal>;
using PointsSectionDF_BlockSparse = PointsSectionDF<csBlockMatrix::BlockSparse>;
using PointsSectionDF_BlockCoeffSparse = PointsSectionDF<csBlockMatrix::BlockCoeffSparse>;

//-----------------------------------------------------------------------------------------------------------------

template <int matType>
class PointsSectionDD : public csNelson::SingleSection<PointsSectionDD<matType>, Point, matType, double, csBlockMatrix::Dynamic,
                                                       csBlockMatrix::Dynamic> {
  std::vector<Point> _points;
  Point _fixedPoint;

public:
  PointsSectionDD() {
    _points.resize(numBlocks);
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  int parameterSize() const override {
    return secSizeFix;
  }

  int numParameters() const override {
    return _points.size();
  }
  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionDD_BlockDense = PointsSectionDD<csBlockMatrix::BlockDense>;
using PointsSectionDD_BlockDiagonal = PointsSectionDD<csBlockMatrix::BlockDiagonal>;
using PointsSectionDD_BlockSparse = PointsSectionDD<csBlockMatrix::BlockSparse>;
using PointsSectionDD_BlockCoeffSparse = PointsSectionDD<csBlockMatrix::BlockCoeffSparse>;

//-----------------------------------------------------------------------------------------------------------------

template <int matType>
class PointsSectionVF
    : public csNelson::SingleSection<PointsSectionVF<matType>, Point, matType, double, csBlockMatrix::Variable, numBlocks> {
  std::array<Point, numBlocks> _points;
  std::vector<int> _sizes;
  Point _fixedPoint;

public:
  PointsSectionVF()
      : _sizes(numBlocks, secSizeFix) {
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  const std::vector<int>& parameterSize() const override {
    return _sizes;
  }
  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionVF_BlockDense = PointsSectionVF<csBlockMatrix::BlockDense>;
using PointsSectionVF_BlockDiagonal = PointsSectionVF<csBlockMatrix::BlockDiagonal>;
using PointsSectionVF_BlockSparse = PointsSectionVF<csBlockMatrix::BlockSparse>;
using PointsSectionVF_BlockCoeffSparse = PointsSectionVF<csBlockMatrix::BlockCoeffSparse>;

//-----------------------------------------------------------------------------------------------------------------

template <int matType>
class PointsSectionVD : public csNelson::SingleSection<PointsSectionVD<matType>, Point, matType, double, csBlockMatrix::Variable,
                                                       csBlockMatrix::Dynamic> {
  std::vector<Point> _points;
  std::vector<int> _sizes;
  Point _fixedPoint;

public:
  PointsSectionVD()
      : _sizes(numBlocks, secSizeFix) {
    _points.resize(numBlocks);
    this->parametersReady();
  }

  virtual const Point& parameter(csNelson::NodeId i) const {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }
  virtual Point& parameter(csNelson::NodeId i) {
    if (i.isVariable())
      return _points[i.id()];
    else {
      assert(i.id() == 0);
      return _fixedPoint;
    }
  }

  const std::vector<int>& parameterSize() const override {
    return _sizes;
  }
  int numFixedParameters() const override {
    return 1;
  }
};

using PointsSectionVD_BlockDense = PointsSectionVD<csBlockMatrix::BlockDense>;
using PointsSectionVD_BlockDiagonal = PointsSectionVD<csBlockMatrix::BlockDiagonal>;
using PointsSectionVD_BlockSparse = PointsSectionVD<csBlockMatrix::BlockSparse>;
using PointsSectionVD_BlockCoeffSparse = PointsSectionVD<csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("SingleSection", "[SingleSection]", PointsSectionFF_BlockDense, PointsSectionFF_BlockDiagonal,
                   PointsSectionFF_BlockSparse, PointsSectionFF_BlockCoeffSparse, PointsSectionFD_BlockDense,
                   PointsSectionFD_BlockDiagonal, PointsSectionFD_BlockSparse, PointsSectionFD_BlockCoeffSparse,
                   PointsSectionDF_BlockDense, PointsSectionDF_BlockDiagonal, PointsSectionDF_BlockSparse,
                   PointsSectionDF_BlockCoeffSparse, PointsSectionDD_BlockDense, PointsSectionDD_BlockDiagonal,
                   PointsSectionDD_BlockSparse, PointsSectionDD_BlockCoeffSparse, PointsSectionVF_BlockDense,
                   PointsSectionVF_BlockDiagonal, PointsSectionVF_BlockSparse, PointsSectionVF_BlockCoeffSparse,
                   PointsSectionVD_BlockDense, PointsSectionVD_BlockDiagonal, PointsSectionVD_BlockSparse,
                   PointsSectionVD_BlockCoeffSparse) {
  std::cout << "-------------------------------------------------------" << std::endl;
  TestType pss;
  //  REQUIRE(pss.parameterSize() == secSizeFix);
  REQUIRE(pss.numParameters() == numBlocks);

  // unary edge
  for (int i = 0; i < numBlocks; i++) {
    pss.addEdge(i, new EdgeUnaryTest<TestType>(i));
  }

  // unary edge on fixed (not really useful...)
  pss.addEdge(csNelson::NodeId::fixed(0), new EdgeUnaryTest<TestType>(0));

  // add "odometry"  edges (not if mat diagonal)
  if (pss.matType() != csBlockMatrix::BlockDiagonal) {
    for (int i = 0; i < numBlocks - 1; i++) {
      pss.addEdge(i, i + 1, new EdgeBinaryTest<TestType>(i, i + 1));
    }
  }

  // add fixed odometry edge
  pss.addEdge(csNelson::NodeId::fixed(0), 0, new EdgeBinaryTest<TestType>(0, 0));
  pss.addEdge(numBlocks - 1, csNelson::NodeId::fixed(0), new EdgeBinaryTest<TestType>(numBlocks - 1, 0));

  // add "extreme" bin relation (not if mat diagonal)
  if (pss.matType() != csBlockMatrix::BlockDiagonal) {
    for (int i = 0; i < numBlocks - 1; i++) {
      pss.addEdge(i, numBlocks - 1, new EdgeBinaryTest<TestType>(i, numBlocks - 1));
    }
  }

  // ready!
  pss.structureReady();

  for (int i = 0; i < numBlocks; i++) {
    REQUIRE(pss.sparsityPattern().has(i, i));
  }
  // check "odometry"  edges (not if mat diagonal)
  if (pss.matType() != csBlockMatrix::BlockDiagonal) {
    for (int i = 0; i < numBlocks; i++) {
      REQUIRE(pss.sparsityPattern().has(i, i));
    }
    for (int i = 0; i < numBlocks - 1; i++) {
      REQUIRE(pss.sparsityPattern().has(i, i + 1));
    }
  }

  // check "extreme" bin relation (not if mat diagonal)
  if (pss.matType() != csBlockMatrix::BlockDiagonal) {
    for (int i = 0; i < numBlocks - 1; i++) {
      REQUIRE(pss.sparsityPattern().has(i, numBlocks - 1));
    }
  }

  pss.update(true);
}
