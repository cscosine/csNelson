#pragma once
#include <csNelson/Global.hpp>

#include <csBlockMatrix/Global.hpp>
#include <csBlockMatrix/MatrixTypeTraits.hpp>
#include <csBlockMatrix/SparsityPattern.hpp>
#include <csBlockMatrix/VectorBlock.hpp>

namespace csNelson {

/*
U  W
W' V
*/

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv = csBlockMatrix::Dynamic,
          int NBVv = csBlockMatrix::Dynamic>
class DoubleSectionHessianMatrices {
public:
  using MatTraitsU = csBlockMatrix::MatrixBlockIterableTypeTraits<matTypeUv, Tv, csBlockMatrix::ColMajor, BUv, BUv, NBUv, NBUv>;
  using MatTraitsV = csBlockMatrix::MatrixBlockIterableTypeTraits<matTypeVv, Tv, csBlockMatrix::ColMajor, BVv, BVv, NBVv, NBVv>;
  using MatTraitsW = csBlockMatrix::MatrixBlockIterableTypeTraits<matTypeWv, Tv, csBlockMatrix::RowMajor, BUv, BVv, NBUv, NBVv>;

  using MatTypeU = typename MatTraitsU::MatrixType;
  using MatTypeV = typename MatTraitsV::MatrixType;
  using MatTypeW = typename MatTraitsW::MatrixType;

private:
  MatTypeU _U;
  MatTypeV _V;
  MatTypeW _W;

public:
  DoubleSectionHessianMatrices();
  virtual ~DoubleSectionHessianMatrices();
  inline MatTypeU& U() {
    return _U;
  }
  inline MatTypeV& V() {
    return _V;
  }
  inline MatTypeW& W() {
    return _W;
  }
  inline const MatTypeU& U() const {
    return _U;
  }
  inline const MatTypeV& V() const {
    return _V;
  }
  inline const MatTypeW& W() const {
    return _W;
  }
};

template <class Tv, int BUv, int BVv, int NBUv = csBlockMatrix::Dynamic, int NBVv = csBlockMatrix::Dynamic>
class DoubleSectionHessianVectors {
public:
  using VecTypeU = csBlockMatrix::VectorBlock<Tv, BUv, NBUv>;
  using VecTypeV = csBlockMatrix::VectorBlock<Tv, BVv, NBVv>;

private:
  VecTypeU _bU;
  VecTypeV _bV;

public:
  DoubleSectionHessianVectors();
  virtual ~DoubleSectionHessianVectors();
  inline VecTypeU& bU() {
    return _bU;
  }
  inline VecTypeV& bV() {
    return _bV;
  }
  inline const VecTypeU& bU() const {
    return _bU;
  }
  inline const VecTypeV& bV() const {
    return _bV;
  }

  void setZero() {
    _bU.setZero();
    _bV.setZero();
  }
};

} // namespace csNelson
