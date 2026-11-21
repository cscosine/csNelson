#pragma once
#include <csNelson/Global.hpp>

#include <csBlockMatrix/Global.hpp>
#include <csBlockMatrix/MatrixTypeTraits.hpp>
#include <csBlockMatrix/SparsityPattern.hpp>
#include <csBlockMatrix/VectorBlock.hpp>

namespace csNelson {

template <int matTypeV, class Tv, int Bv, int NBv> struct SingleSectionHessianTraits {
  static constexpr int matType = matTypeV;
  using Type = Tv;
  static constexpr int B = Bv;
  static constexpr int NB = NBv;
};

template <int matTypeV, class T, int B, int NB = csBlockMatrix::Dynamic> class SingleSectionHessian {
public:
  using MatTraits = csBlockMatrix::MatrixBlockIterableTypeTraits<matTypeV, T, csBlockMatrix::ColMajor, B, B, NB, NB>;

  using Traits = SingleSectionHessianTraits<matTypeV, T, B, NB>;

  using MatType = typename MatTraits::MatrixType;
  using VecType = csBlockMatrix::VectorBlock<T, B, NB>;

  using BlockSizeTypePar = typename MatType::RowTraits::BlockSizeTypePar;

private:
  MatType _H;
  VecType _b;
  double _chi2;

public:
  SingleSectionHessian();
  virtual ~SingleSectionHessian();

  void resize(BlockSizeTypePar blockSizes, int nBlocks,
              const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &sp);

  inline void setChi2(double chi2) { this->_chi2 = chi2; }

  inline double chi2() const { return _chi2; }

  inline const MatType &H() const { return _H; }

  inline MatType &H() { return _H; }

  inline const VecType &b() const { return _b; }

  inline VecType &b() { return _b; }

  void clearAll();

  T maxAbsValBVect() const;
};
} // namespace csNelson
