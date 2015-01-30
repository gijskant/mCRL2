// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/pbes2data_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_DETAIL_PBES2DATA_REWRITER_H
#define MCRL2_PBES_DETAIL_PBES2DATA_REWRITER_H

#include "mcrl2/pbes/builder.h"

namespace mcrl2 {

namespace pbes_system {

namespace detail {

// transforms outer level pbes operators to their data equivalents, for the following operators:
// x && y
// x || y
// x => y
// exists d:D. x
// forall d:D. x
template <typename Derived>
struct pbes2data_builder: public pbes_expression_builder<Derived>
{
  typedef pbes_expression_builder<Derived> super;
  using super::enter;
  using super::leave;
  using super::apply;

  bool is_data(const pbes_expression& x) const
  {
    return pbes_system::is_data(x);
  }

  bool is_not(const pbes_expression& x) const
  {
    return pbes_system::is_not(x);
  }

  bool is_and(const pbes_expression& x) const
  {
    return pbes_system::is_and(x);
  }

  bool is_or(const pbes_expression& x) const
  {
    return pbes_system::is_or(x);
  }

  bool is_imp(const pbes_expression& x) const
  {
    return pbes_system::is_imp(x);
  }

  bool is_forall(const pbes_expression& x) const
  {
    return pbes_system::is_forall(x);
  }

  bool is_exists(const pbes_expression& x) const
  {
    return pbes_system::is_exists(x);
  }

  pbes_expression arg(const pbes_expression& x) const
  {
    return accessors::arg(x);
  }

  pbes_expression left(const pbes_expression& x) const
  {
    return accessors::left(x);
  }

  pbes_expression right(const pbes_expression& x) const
  {
    return accessors::right(x);
  }

  // transforms outer level data operators to their pbes equivalents, for the following operators:
  // !x
  // x && y
  // x || y
  // x => y
  // exists d:D. x
  // forall d:D. x
  pbes_expression pbes2data(const pbes_expression& x) const
  {
    if (is_not(x))
    {
      pbes_expression operand = pbes2data(arg(x));
      if (is_data(operand))
      {
        return data::sort_bool::not_(data::data_expression(operand));
      }
      else
      {
        return pbes_system::optimized_not(operand);
      }
    }
    else if (is_and(x))
    {
      pbes_expression l = pbes2data(left(x));
      pbes_expression r = pbes2data(right(x));
      if (is_data(l) && is_data(r))
      {
        return data::sort_bool::and_(data::data_expression(l), data::data_expression(r));
      }
      else
      {
        return pbes_system::optimized_and(l, r);
      }
    }
    else if (is_or(x))
    {
      pbes_expression l = pbes2data(left(x));
      pbes_expression r = pbes2data(right(x));
      if (is_data(l) && is_data(r))
      {
        return data::sort_bool::or_(data::data_expression(l), data::data_expression(r));
      }
      else
      {
        return pbes_system::optimized_or(l, r);
      }

    }
    else if (is_imp(x))
    {
      pbes_expression l = pbes2data(left(x));
      pbes_expression r = pbes2data(right(x));
      if (is_data(l) && is_data(r))
      {
        return data::sort_bool::implies(data::data_expression(l), data::data_expression(r));
      }
      else
      {
        return pbes_system::optimized_imp(l, r);
      }

    }
    else if (is_forall(x))
    {
      pbes_system::forall y(x);
      pbes_expression operand = pbes2data(arg(y));
      if (is_data(operand))
      {
        return data::forall(y.variables(), data::data_expression(operand));
      }
      else
      {
        return pbes_system::optimized_forall(y.variables(), operand);
      }
    }
    else if (is_exists(x))
    {
      pbes_system::exists y(x);
      pbes_expression operand = pbes2data(arg(y));
      if (is_data(operand))
      {
        return data::exists(y.variables(), data::data_expression(operand));
      }
      else
      {
        return pbes_system::optimized_exists(y.variables(), operand);
      }
    }
    return x;
  }

  pbes_expression apply(const pbes_expression& x)
  {
    return pbes2data(x);
  }
};

template <typename T>
T pbes2data(const T& x,
            typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0
           )
{
  return core::make_apply_builder<pbes2data_builder>().apply(x);
}

template <typename T>
void pbes2data(T& x,
               typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0
              )
{
  core::make_apply_builder<pbes2data_builder>().apply(x);
}

} // namespace detail

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_DETAIL_PBES2DATA_REWRITER_H
