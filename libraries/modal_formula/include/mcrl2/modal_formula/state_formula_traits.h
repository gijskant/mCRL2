// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/state_formula_traits.h
/// \brief Contains term traits for state_formula.

#ifndef MCRL2_STATE_FORMULA_TERM_TRAITS_H
#define MCRL2_STATE_FORMULA_TERM_TRAITS_H

#include <functional>
#include "mcrl2/core/term_traits.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/data/variable.h"
#include "mcrl2/data/detail/data_sequence_algorithm.h"
#include "mcrl2/utilities/exception.h"
#include "mcrl2/utilities/detail/join.h"

namespace mcrl2
{

namespace core
{

/// \brief Contains type information for data expressions.
template <>
struct term_traits<state_formulas::state_formula>
{
  /// \brief The term type
  typedef state_formulas::state_formula term_type;

  /// \brief The variable type
  typedef data::variable variable_type;

  /// \brief The variable sequence type
  typedef data::variable_list variable_sequence_type;

  /// \brief The value true
  /// \return The value true
  static inline
  term_type true_()
  {
    return state_formulas::true_();
  }

  /// \brief The value false
  /// \return The value false
  static inline
  term_type false_()
  {
    return state_formulas::false_();
  }

  /// \brief Operator not
  /// \param p A term
  /// \return Operator not applied to p
  static inline
  term_type not_(const term_type &p)
  {
    return state_formulas::not_(p);
  }

  /// \brief Operator and
  /// \param p A term
  /// \param q A term
  /// \return Operator and applied to p and q
  static inline
  term_type and_(const term_type& p, const term_type& q)
  {
    return state_formulas::and_(p, q);
  }

  /// \brief Operator or
  /// \param p A term
  /// \param q A term
  /// \return Operator or applied to p and q
  static inline
  term_type or_(const term_type& p, const term_type& q)
  {
    return state_formulas::or_(p, q);
  }

  /// \brief Operator imp
  /// \param p A term
  /// \param q A term
  /// \return Operator or applied to p and q
  static inline
  term_type imp(const term_type& p, const term_type& q)
  {
    return state_formulas::imp(p, q);
  }

  /// \brief Operator forall
  /// \param d A sequence of variables
  /// \param p A term
  /// \return Operator forall applied to d and p
  static inline
  term_type forall(const variable_sequence_type &d, const term_type& p)
  {
    return state_formulas::forall(d, p);
  }

  /// \brief Operator exists
  /// \param d A sequence of variables
  /// \param p A term
  /// \return Operator exists applied to d and p
  static inline
  term_type exists(const variable_sequence_type &d, const term_type& p)
  {
    return state_formulas::exists(d, p);
  }

  /// \brief Test for data term
  /// \param t A term
  /// \return True if the term is a data expression
  static inline
  bool is_data(const term_type& t)
  {
    return data::is_data_expression(t);
  }

  /// \brief Test for value true
  /// \param t A term
  /// \return True if the term has the value true
  static inline
  bool is_true(const term_type& t)
  {
    return state_formulas::is_true(t);
  }

  /// \brief Test for value false
  /// \param t A term
  /// \return True if the term has the value false
  static inline
  bool is_false(const term_type& t)
  {
    return state_formulas::is_false(t);
  }

  /// \brief Test for operator not
  /// \param t A term
  /// \return True if the term is of type not
  static inline
  bool is_not(const term_type& t)
  {
    return state_formulas::is_not(t);
  }

  /// \brief Test for operator and
  /// \param t A term
  /// \return True if the term is of type and
  static inline
  bool is_and(const term_type& t)
  {
    return state_formulas::is_and(t);
  }

  /// \brief Test for operator or
  /// \param t A term
  /// \return True if the term is of type or
  static inline
  bool is_or(const term_type& t)
  {
    return state_formulas::is_or(t);
  }

  /// \brief Test for implication
  /// \param t A term
  /// \return True if the term is an implication
  static inline
  bool is_imp(const term_type& t)
  {
    return state_formulas::is_imp(t);
  }

  /// \brief Test for universal quantification
  /// \param t A term
  /// \return True if the term is an universal quantification
  static inline
  bool is_forall(const term_type& t)
  {
    return state_formulas::is_forall(t);
  }

  /// \brief Test for existential quantification
  /// \param t A term
  /// \return True if the term is an existential quantification
  static inline
  bool is_exists(const term_type& t)
  {
    return state_formulas::is_exists(t);
  }

  /// \brief Test if a term is a variable
  /// \param t A term
  /// \return True if the term is a variable
  static inline
  bool is_variable(const term_type& t)
  {
    return state_formulas::is_variable(t);
  }

  static inline
  term_type left(const term_type& t)
  {
    return atermpp::down_cast<const term_type>(t[0]);
  }

  static inline
  term_type right(const term_type& t)
  {
    return atermpp::down_cast<const term_type>(t[1]);
  }

  static inline
  term_type not_arg(const term_type& t)
  {
    return atermpp::down_cast<const term_type>(t[0]);
  }

  /// \brief Pretty print function
  /// \param t A term
  /// \return A pretty print representation of the term
  static inline
  std::string pp(const term_type& t)
  {
    return state_formulas::pp(t);
  }
};

} // namespace core

} // namespace mcrl2

#endif // MCRL2_STATE_FORMULA_TERM_TRAITS_H
