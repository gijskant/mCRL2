/*
 * modal_formula_optimized.h
 *
 *  Created on: 6 Sep 2013
 *      Author: kant
 */

#ifndef MODAL_FORMULA_OPTIMIZED_H_
#define MODAL_FORMULA_OPTIMIZED_H_

#include <iostream> // for debugging

#include <string>
#include <cassert>
#include "mcrl2/atermpp/aterm_appl.h"
#include "mcrl2/core/print.h"
#include "mcrl2/core/term_traits.h"
#include "mcrl2/core/detail/precedence.h"
//#include "mcrl2/core/detail/struct_core.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/utilities/detail/join.h"
#include "mcrl2/data/optimized_boolean_operators.h"

namespace mcrl2 {

namespace core {

/// \brief Contains type information for state formulas.
template <>
struct term_traits_optimized<state_formulas::state_formula>: public core::term_traits<state_formulas::state_formula>
{
  typedef core::term_traits<state_formulas::state_formula> super;

  /// \brief Make a negation
  /// \param p A state formula
  /// \return The value <tt>!p</tt>
  static inline
  term_type not_(const term_type& p)
  {
    return data::optimized_not(p);
  }

  /// \brief Make a conjunction
  /// \param p A state formula
  /// \param q A state formula
  /// \return The value <tt>p && q</tt>
  static inline
  term_type and_(const term_type& p, const term_type& q)
  {
    return data::optimized_and(p, q);
  }

  /// \brief Make a disjunction
  /// \param p A state formula
  /// \param q A state formula
  /// \return The value <tt>p || q</tt>
  static inline
  term_type or_(const term_type& p, const term_type& q)
  {
    return data::optimized_or(p, q);
  }

  /// \brief Make an implication
  /// \param p A state formula
  /// \param q A state formula
  /// \return The value <tt>p => q</tt>
  static inline
  term_type imp(const term_type& p, const term_type& q)
  {
    return data::optimized_imp(p, q);
  }

  /// \brief Returns or applied to the sequence of state formulas [first, last)
  /// \param first Start of a sequence of state formulas
  /// \param last End of a sequence of state formulas
  /// \return Or applied to the sequence of state formulas [first, last)
  template <typename FwdIt>
  static inline
  term_type join_or(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, or_, false_());
  }

  /// \brief Returns and applied to the sequence of state formulas [first, last)
  /// \param first Start of a sequence of state formulas
  /// \param last End of a sequence of state formulas
  /// \return And applied to the sequence of state formulas [first, last)
  template <typename FwdIt>
  static inline
  term_type join_and(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, and_, true_());
  }

  /// \brief Make a universal quantification
  /// If l is empty, p is returned.
  /// \param l A sequence of data variables
  /// \param p A state formula
  /// \return The value <tt>forall l.p</tt>
  static inline
  term_type forall(const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      return p;
    }
    if (is_false(p))
    {
      // N.B. Here we use the fact that mCRL2 data types are never empty.
      return false_();
    }
    if (is_true(p))
    {
      return true_();
    }
    return state_formulas::forall(l, p);
  }

  /// \brief Make an existential quantification
  /// If l is empty, p is returned.
  /// \param l A sequence of data variables
  /// \param p A state formulas
  /// \return The value <tt>exists l.p</tt>
  static inline
  term_type exists(const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      return p;
    }
    if (is_false(p))
    {
      return false_();
    }
    if (is_true(p))
    {
      // N.B. Here we use the fact that mCRL2 data types are never empty.
      return true_();
    }
    return state_formulas::exists(l, p);
  }

  /// \brief Make an existential modal formula (may).
  /// If formula or operand is false, false is returned.
  /// \param formula A regular formula or action formula.
  /// \param operand The state formula after the action.
  /// \return The value <tt><formula>operand</tt>
  static inline
  term_type may(const regular_formulas::regular_formula& formula, const term_type& operand)
  {
    if (is_false(operand))
    {
      return false_();
    }
    if (action_formulas::is_false(formula))
    {
      return false_();
    }
    return state_formulas::may(formula, operand);
  }

  /// \brief Make an universal modal formula (must).
    /// If formula is false or operand is true, true is returned
    /// \param formula A regular formula or action formula.
    /// \param operand The state formula after the action.
    /// \return The value <tt><formula>operand</tt>
    static inline
    term_type must(const regular_formulas::regular_formula& formula, const term_type& operand)
    {
      if (is_true(operand))
      {
        return true_();
      }
      if (action_formulas::is_false(formula))
      {
        return true_();
      }
      return state_formulas::must(formula, operand);
    }

};

} // namespace core

} // namespace mcrl2

#endif /* MODAL_FORMULA_OPTIMIZED_H_ */
