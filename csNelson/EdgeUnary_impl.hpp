#pragma once
#include <csNelson/EdgeSectionBase_impl.hpp>
#include <csNelson/EdgeUnary.hpp>

namespace csNelson {

template <class Section>
EdgeUnarySectionBase<Section>::EdgeUnarySectionBase() : EdgeUnaryBase(), EdgeSectionBase<Section>() {}

template <class Section> EdgeUnarySectionBase<Section>::~EdgeUnarySectionBase() {}

} // namespace csNelson
