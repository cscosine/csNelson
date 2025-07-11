#pragma once
#include "Global.h"

#include "csBlockMatrix/MatrixTypeTraits.h"
#include "csBlockMatrix/VectorBlock.h"

#include "ParallelExecHelper.h"

namespace csNelson {

struct MatrixDiagInvSettings {
  ParallelExecSettings blockInversion;
  ParallelExecSettings rightVectorMult;
};

template <int matType, class T, int BV, int NBV, int matTypeOut>
class MatrixDiagInv {
public:
  using MatTraitsV = csBlockMatrix::MatrixBlockIterableTypeTraits<matType, T, csBlockMatrix::ColMajor, BV, BV, NBV, NBV>;
  using MatTypeV = typename MatTraitsV::MatrixType;

  using MatTraitsVinv = csBlockMatrix::MatrixBlockIterableTypeTraits<matTypeOut, T, csBlockMatrix::ColMajor, BV, BV, NBV, NBV>;
  using MatTypeVinv = typename MatTraitsVinv::MatrixType;

private:
  MatTypeVinv _Vinv;

  MatrixDiagInvSettings _settings;

public:
  const MatrixDiagInvSettings& settings() const {
    return _settings;
  }
  MatrixDiagInvSettings& settings() {
    return _settings;
  }

  void init(const MatTypeV& V);
  void compute(const MatTypeV& V, T relLambda, T absLambda);

  MatTypeVinv& Vinv() {
    return _Vinv;
  }

  void rightMultVector(const csBlockMatrix::VectorBlock<T, BV, NBV>& v, csBlockMatrix::VectorBlock<T, BV, NBV>& res) const;
};

} // namespace csNelson
