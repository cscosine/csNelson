#pragma once
#include "EdgeBinary.h"

namespace csNelson {

template <class Section>
EdgeBinarySectionBase<Section>::EdgeBinarySectionBase()
    : EdgeBinaryBase()
    , EdgeSectionBase<Section>() {}

template <class Section>
EdgeBinarySectionBase<Section>::~EdgeBinarySectionBase() {}

} // namespace csNelson
