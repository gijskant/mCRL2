// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/gauss_elimination.h
/// \brief Gauss elimination algorithm for pbes equation systems.

#ifndef MCRL2_PBES_GAUSS_ELIMINATION_H
#define MCRL2_PBES_GAUSS_ELIMINATION_H

#include <sstream>

#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/replace.h"
#include "mcrl2/pbes/rewriter.h"
#include "mcrl2/utilities/logger.h"

namespace mcrl2
{

namespace pbes_system
{

/// \brief Algorithm class for the Gauss elimination algorithm for solving
/// systems of (P)BES equations.

template <typename ExpressionTraits>
class gauss_elimination_algorithm
{
  public:
    typedef typename ExpressionTraits::expression_type expression_type;
    typedef typename ExpressionTraits::variable_type variable_type;
    typedef typename ExpressionTraits::equation_type equation_type;

  protected:

    std::string print_equation(const equation_type& eq) const
    {
      return ExpressionTraits::print(eq);
    }

    template <typename Iter>
    std::string print_equations(Iter first, Iter last) const
    {
      std::ostringstream out;
      for (Iter i = first; i != last; ++i)
      {
        out << "  " << print_equation(*i) << std::endl;
      }
      return out.str();
    }

  public:
    /// \brief Runs the algorithm. Applies Gauss elimination to the sequence of pbes equations [first, last).
    /// \param first Start of a range of pbes equations
    /// \param last End of a range of pbes equations
    /// \param solve An equation solver

    template <typename Iter, typename FixpointEquationSolver>
    void run(Iter first, Iter last, FixpointEquationSolver solve)
    {
      mCRL2log(log::debug) << "equations before solving\n" << print_equations(first, last);
      if (first == last)
      {
        return;
      }

      Iter i = last;
      while (i != first)
      {
        --i;
        mCRL2log(log::verbose) << " * solving equation " << i->variable().name() << std::endl;
        mCRL2log(log::debug) << "solving equation\n  before: " << print_equation(*i);
        solve(*i);
        mCRL2log(log::debug) << "   after: " << print_equation(*i) << "\n";
        for (Iter j = first; j != i; ++j)
        {
          j->formula() = ExpressionTraits::substitute(j->formula(), i->variable(), i->formula());
        }
        mCRL2log(log::debug1) << "equations after substitution\n" << print_equations(first, last);
      }
      mCRL2log(log::debug) << "equations after solving\n" << print_equations(first, last);
    }
};

/// \brief Approximation algorithm

template <typename BooleanExpressionTraits, typename Compare, typename Rewr>
struct approximate
{
  Compare m_compare;
  size_t m_max_iterations;
  Rewr m_rewr;

  approximate(Compare compare, size_t max_iterations = 0, Rewr rewr = 0)
    : m_compare(compare),
      m_max_iterations(max_iterations),
      m_rewr(rewr)
  {
  }

  void operator()(typename BooleanExpressionTraits::equation_type& eq) const
  {
    typedef BooleanExpressionTraits tr;
    typedef typename BooleanExpressionTraits::expression_type expression_type;
    typedef typename BooleanExpressionTraits::variable_type variable_type;

    const expression_type& phi = eq.formula();
    const variable_type& X = eq.variable();
    expression_type next = eq.symbol() == tr::nu() ? tr::true_() : tr::false_();
    expression_type prev;
    size_t i = 0;
    do
    {
      mCRL2log(log::verbose) << "   - iteration " << i << " / " << m_max_iterations << std::endl;
      prev = next;
      next = tr::substitute(phi, X, prev);
      if (m_rewr != 0) { next = (*m_rewr)(next); }
      i++;
    }
    while (!m_compare(prev, next) && (m_max_iterations==0 || i < m_max_iterations));
    eq.formula() = prev;
  }
};

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_GAUSS_ELIMINATION_H
