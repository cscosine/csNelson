#pragma once
#include <csNelson/Global.hpp>

#include <csBlockMatrix/Global.hpp>
#include <csBlockMatrix/MatrixTypeTraits.hpp>
#include <csBlockMatrix/SparsityPattern.hpp>
#include <csBlockMatrix/VectorBlock.hpp>

#include <csNelson/DoubleSectionHessianMatrices.hpp>

namespace csNelson {

/*
U  W
W' V
*/

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
struct DoubleSectionHessianTraits {
  static constexpr int matTypeU = matTypeUv;
  static constexpr int matTypeW = matTypeWv;
  static constexpr int matTypeV = matTypeVv;

  using Type = Tv;

  static constexpr int BU = BUv;
  static constexpr int BV = BVv;
  static constexpr int NBU = NBUv;
  static constexpr int NBV = NBVv;

  using HessianMatricesType = DoubleSectionHessianMatrices<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv, NBVv>;
  using HessianVectorsType = DoubleSectionHessianVectors<Tv, BUv, BVv, NBUv, NBVv>;
};

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv = csBlockMatrix::Dynamic,
          int NBVv = csBlockMatrix::Dynamic>
class DoubleSectionHessian {
public:
  using Traits = DoubleSectionHessianTraits<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv, NBVv>;

  using HessianMatricesT = typename Traits::HessianMatricesType;
  using HessianVectorsT = typename Traits::HessianVectorsType;

  using MatTraitsU = typename HessianMatricesT::MatTraitsU;
  using MatTraitsV = typename HessianMatricesT::MatTraitsV;
  using MatTraitsW = typename HessianMatricesT::MatTraitsW;

  using MatTypeU = typename MatTraitsU::MatrixType;
  using MatTypeV = typename MatTraitsV::MatrixType;
  using MatTypeW = typename MatTraitsW::MatrixType;

  using VecTypeU = typename HessianVectorsT::VecTypeU;
  using VecTypeV = typename HessianVectorsT::VecTypeV;

  using BlockSizeTypeParU = typename MatTypeU::RowTraits::BlockSizeTypePar;
  using BlockSizeTypeParV = typename MatTypeV::RowTraits::BlockSizeTypePar;

private:
  HessianMatricesT _H;

  HessianVectorsT _b;

  double _chi2;

public:
  DoubleSectionHessian();
  virtual ~DoubleSectionHessian();

  void resize(BlockSizeTypeParU blockSizesU, int nBlocksU, BlockSizeTypeParV blockSizesV, int nBlocksV,
              const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &spU,
              const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &spV,
              const csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>::CSPtr &spW);

  HessianMatricesT &H() { return _H; }
  const HessianMatricesT &H() const { return _H; }

  HessianVectorsT &b() { return _b; }
  const HessianVectorsT &b() const { return _b; }

  inline void setChi2(double chi2) { this->_chi2 = chi2; }

  inline double chi2() const { return _chi2; }

  inline const MatTypeU &U() const { return _H.U(); }
  inline const MatTypeV &V() const { return _H.V(); }
  inline const MatTypeW &W() const { return _H.W(); }

  inline MatTypeU &U() { return _H.U(); }
  inline MatTypeV &V() { return _H.V(); }
  inline MatTypeW &W() { return _H.W(); }

  inline const VecTypeU &bU() const { return _b.bU(); }
  inline const VecTypeV &bV() const { return _b.bV(); }

  inline VecTypeU &bU() { return _b.bU(); }
  inline VecTypeV &bV() { return _b.bV(); }

  void clearAll();

  Tv maxAbsValBVect() const;
};

} // namespace csNelson
