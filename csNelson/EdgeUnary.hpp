#pragma once
#include "EdgeSectionBase.hpp"
#include "EdgeUnary.h"

namespace csNelson {

template <class Section>
EdgeUnarySectionBase<Section>::EdgeUnarySectionBase()
    : EdgeUnaryBase()
    , EdgeSectionBase<Section>() {}

template <class Section>
EdgeUnarySectionBase<Section>::~EdgeUnarySectionBase() {}

} // namespace csNelson
