// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/remove_parameters.h
/// \brief Functions for removing insignificant parameters from state formulas.

#ifndef MCRL2_MODAL_FORMULA_REMOVE_PARAMETERS_H
#define MCRL2_MODAL_FORMULA_REMOVE_PARAMETERS_H

#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/state_formula.h"
#include <algorithm>
#include <functional>
#include <map>
#include <vector>

namespace mcrl2
{

namespace state_formulas
{

/// \cond INTERNAL_DOCS
namespace detail
{

/// \brief Removes elements with indices in a given sequence from the sequence l
/// This is a copy from the remove_elements function in mcrl2/pbes/remove_parameters.h.
/// \param l A sequence of terms
/// \param to_be_removed A sequence of integers
/// \return The removal result
template <typename Term>
atermpp::term_list<Term> remove_elements(atermpp::term_list<Term> l, const std::vector<size_t>& to_be_removed)
{
  size_t index = 0;
  std::vector<Term> result;
  std::vector<size_t>::const_iterator j = to_be_removed.begin();
  for (typename atermpp::term_list<Term>::iterator i = l.begin(); i != l.end(); ++i, ++index)
  {
    if (j != to_be_removed.end() && index == *j)
    {
      ++j;
    }
    else
    {
      result.push_back(*i);
    }
  }
  return atermpp::term_list< Term >(result.begin(),result.end());
}

template <typename Derived>
struct remove_parameters_builder: public state_formula_builder<Derived>
{
  typedef state_formula_builder<Derived> super;
  using super::enter;
  using super::leave;
  using super::apply;

  const std::map<core::identifier_string, std::vector<size_t> >& to_be_removed;

  remove_parameters_builder(const std::map<core::identifier_string, std::vector<size_t> >& to_be_removed_)
    : to_be_removed(to_be_removed_)
  {}

  inline std::string print_vector(std::vector<size_t> x)
  {
    std::stringstream s;
    bool first = true;
    for (auto it = x.begin(); it != x.end(); ++it)
    {
      if (!first)
      {
        s << ", ";
      }
      s << (*it);
      first = false;
    }
    return s.str();
  }

  // to prevent default apply being called
  data::data_expression apply(const data::data_expression& x)
  {
  	return x;
  }

  state_formulas::variable apply(const state_formulas::variable& x)
  {
    mCRL2log(log::debug) << "variable " << x.name() << std::endl;
    std::map<core::identifier_string, std::vector<size_t> >::const_iterator i = to_be_removed.find(x.name());
    if (i != to_be_removed.end())
    {
      mCRL2log(log::debug) << " - removing elements: " << print_vector(i->second) << std::endl;
      return state_formulas::variable(x.name(), detail::remove_elements(x.arguments(), i->second));
    }
    return x;
  }

  state_formulas::nu apply(const state_formulas::nu& x)
  {
    mCRL2log(log::debug) << "nu " << x.name() << std::endl;
    state_formulas::state_formula x_operand = static_cast<Derived&>(*this).apply(x.operand());
    std::map<core::identifier_string, std::vector<size_t> >::const_iterator i = to_be_removed.find(x.name());
    if (i == to_be_removed.end())
    {
      return state_formulas::nu(x.name(), x.assignments(), x_operand);
    }
    else
    {
      mCRL2log(log::debug) << " - removing elements: " << print_vector(i->second) << std::endl;
      return state_formulas::nu(x.name(), detail::remove_elements(x.assignments(), i->second), x_operand);
    }
  }

  state_formulas::mu apply(const state_formulas::mu& x)
  {
    mCRL2log(log::debug) << "mu " << x.name() << std::endl;
    state_formulas::state_formula x_operand = static_cast<Derived&>(*this).apply(x.operand());
    std::map<core::identifier_string, std::vector<size_t> >::const_iterator i = to_be_removed.find(x.name());
    if (i == to_be_removed.end())
    {
      return state_formulas::mu(x.name(), x.assignments(), x_operand);
    }
    else
    {
      mCRL2log(log::debug) << " - removing elements: " << print_vector(i->second) << std::endl;
      return state_formulas::mu(x.name(), detail::remove_elements(x.assignments(), i->second), x_operand);
    }
  }

  modal_equation apply(modal_equation& x)
  {
    state_formulas::state_formula x_formula = static_cast<Derived&>(*this).apply(x.formula());
    x.set_formula(x_formula);
    std::map<core::identifier_string, std::vector<size_t> >::const_iterator i = to_be_removed.find(x.name());
    if (i != to_be_removed.end())
    {
      x.set_parameters(detail::remove_elements(x.parameters(), i->second));
    }
    return x;
  }

  modal_equation_system apply(modal_equation_system& x)
  {
    std::vector<modal_equation> equations;
    for(auto e: x.equations())
    {
      (*this).apply(e);
      equations.push_back(e);
    }
    x.set_equations(equations);
    x.set_initial_state((*this).apply(x.initial_state()));
    return x;
  }

};
} // namespace detail
/// \endcond

/// \brief Removes parameters from state formulas
/// \param x A state formula
/// \param to_be_removed A mapping that maps propositional variable names to indices of parameters that are removed
/// \return The expression \p x with parameters removed according to the mapping \p to_be_removed
template <typename T>
void remove_parameters(T& x,
                       const std::map<core::identifier_string, std::vector<size_t> >& to_be_removed)
{
  x = core::make_apply_builder_arg1<detail::remove_parameters_builder>(to_be_removed).apply(x);
}

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_REMOVE_PARAMETERS_H
