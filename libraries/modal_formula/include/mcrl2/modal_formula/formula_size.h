// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/formula_size.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_FORMULA_SIZE_H_
#define MCRL2_MODAL_FORMULAS_FORMULA_SIZE_H_

#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/traverser.h"

namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

// \brief Visitor for computing state formula size.
struct state_formula_size_traverser: public state_formula_traverser<state_formula_size_traverser>
{
  typedef state_formula_traverser<state_formula_size_traverser> super;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  size_t size;

public:
  state_formula_size_traverser():
    size(0)
  {
  }

  size_t get_size()
  {
    return size;
  }

  void enter(const state_formulas::state_formula& x)
  {
    size++;
  }
};

/// \brief Computes the size of the state formula, i.e., the number of subformulas.
size_t formula_size(const state_formulas::state_formula& x)
{
  state_formula_size_traverser f;
  f.apply(x);
  return f.get_size();
}

/// \brief Computes the combined size of the state formulas of the modal equation system, i.e., the number of subformulas.
size_t formula_size(const modal_equation_system& x)
{
  state_formula_size_traverser f;
  for (auto e: x.equations())
  {
    f.apply(e.formula());
  }
  return x.equations().size() + f.get_size();
}

/// \endcond

} // namespace state_formulas

} // namespace mcrl2


#endif /* MCRL2_MODAL_FORMULAS_FORMULA_SIZE_H_ */
