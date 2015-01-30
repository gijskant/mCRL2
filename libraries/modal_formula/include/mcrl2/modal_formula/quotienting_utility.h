// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/quotienting_utils.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_QUOTIENTING_UTILITY_H_
#define MCRL2_MODAL_FORMULAS_QUOTIENTING_UTILITY_H_

namespace mcrl2 {

namespace state_formulas {

struct quotienting_options
{
  quotienting_options() :
    use_rewrite_cache(true),
    unfold_unguarded_recursion(false),
    use_vector_map(false)
  {  }

  bool use_rewrite_cache;
  bool unfold_unguarded_recursion;
  bool use_vector_map;
};

} // namespace state_formulas

} // namespace mcrl2

#endif /* MCRL2_MODAL_FORMULAS_QUOTIENTING_UTILITY_H_ */

