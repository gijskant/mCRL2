// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/transform_modal_equation_system.h
/// \brief Transform state formulas to and from modal equation systems.

#ifndef MCRL2_MODAL_TRANSFORM_MODAL_EQUATION_SYSTEM_H_
#define MCRL2_MODAL_TRANSFORM_MODAL_EQUATION_SYSTEM_H_

#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include "boost/static_assert.hpp"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/modal_formula/state_formulas_optimized.h"
#include "mcrl2/modal_formula/traverser.h"
#include "mcrl2/utilities/detail/join.h"


namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

state_formula transform(const modal_equation_system& x);

modal_equation_system transform(const state_formula& x);

// \brief Visitor for transforming a modal equation system to a state formula.
struct modal_equation_transform_builder: public state_formula_builder<modal_equation_transform_builder>
{
  typedef state_formula_builder<modal_equation_transform_builder> super;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  modal_equation_system m_mes;
  std::stack<std::set<core::identifier_string> > declared_variables_stack;

public:

  modal_equation_transform_builder(const modal_equation_system& mes):
    m_mes(mes)
  {
    declared_variables_stack.push(std::set<core::identifier_string>());
  }

  state_formula apply(const variable& x)
  {
    state_formula result;
    std::set<core::identifier_string> declared_variables = declared_variables_stack.top();
    if (declared_variables.find(x.name()) != declared_variables.end())
    {
      result = x;
    }
    else
    {
      declared_variables.insert(x.name());
      declared_variables_stack.push(declared_variables);
      // expand variable
      modal_equation eq = m_mes.equation_map()[x.name()];
      std::vector<data::assignment> assignments;
      auto it1 = eq.parameters().begin();
      auto it2 = x.arguments().begin();
      for(; it1 != eq.parameters().end() && it2 != x.arguments().end(); ++it1, ++it2)
      {
        assignments.push_back(data::assignment((*it1), (*it2)));
      }
      state_formula operand = (*this).apply(eq.formula());
      if (eq.symbol().is_nu())
      {
        result = state_formulas::nu(x.name(), data::assignment_list(assignments.begin(), assignments.end()), operand);
      }
      else
      {
        result = state_formulas::mu(x.name(), data::assignment_list(assignments.begin(), assignments.end()), operand);
      }
      declared_variables_stack.pop();
    }
    return result;
  }

  state_formula apply(const nu& x)
  {
    throw std::runtime_error("No 'nu'-operator allowed in right hand sides of modal equations.");
  }

  state_formula apply(const mu& x)
  {
    throw std::runtime_error("No 'mu'-operator allowed in right hand sides of modal equations.");
  }

};

// \brief Visitor for transforming a state formula to a modal equation system.
struct state_formula_transform_builder: public state_formula_builder<state_formula_transform_builder>
{
  typedef state_formula_builder<state_formula_transform_builder> super;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  std::vector<modal_equation> m_equations;
  std::stack<std::vector<data::variable>> quantifier_variable_stack;
  std::map<core::identifier_string, std::vector<data::variable>> variable_to_added_arguments_map;

public:

  state_formula_transform_builder()
  {
    quantifier_variable_stack.push(std::vector<data::variable>());
  }

  std::vector<modal_equation> equations()
  {
    return m_equations;
  }

  state_formula apply(const variable& x)
  {
    (*this).enter(x);
    std::vector<data::data_expression> arguments(x.arguments().begin(), x.arguments().end());
    for(auto a : variable_to_added_arguments_map[x.name()])
    {
      arguments.push_back(a);
    }
    variable result(x.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    (*this).leave(x);
    return result;
  }

  state_formula apply(const mu& x)
  {
    (*this).enter(x);
    variable_to_added_arguments_map[x.name()] = quantifier_variable_stack.top();
    state_formula rhs = (*this).apply(x.operand());
    std::vector<data::variable> parameters;
    std::vector<data::data_expression> arguments;
    for(auto a : x.assignments())
    {
      parameters.push_back(a.lhs());
      arguments.push_back(a.rhs());
    }
    for(auto q : quantifier_variable_stack.top())
    {
      parameters.push_back(q);
      arguments.push_back(q);
    }
    modal_equation eq(fixpoint_symbol::mu(), x.name(), data::variable_list(parameters.begin(), parameters.end()), rhs);
    m_equations.push_back(eq);
    variable v(x.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    (*this).leave(x);
    return v;
  }

  state_formula apply(const nu& x)
  {
    (*this).enter(x);
    variable_to_added_arguments_map[x.name()] = quantifier_variable_stack.top();
    state_formula rhs = (*this).apply(x.operand());
    std::vector<data::variable> parameters;
    std::vector<data::data_expression> arguments;
    for(auto a : x.assignments())
    {
      parameters.push_back(a.lhs());
      arguments.push_back(a.rhs());
    }
    for(auto q : quantifier_variable_stack.top())
    {
      parameters.push_back(q);
      arguments.push_back(q);
    }
    modal_equation eq(fixpoint_symbol::nu(), x.name(), data::variable_list(parameters.begin(), parameters.end()), rhs);
    m_equations.push_back(eq);
    variable v(x.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    (*this).leave(x);
    return v;
  }

  state_formula apply(const forall& x)
  {
    (*this).enter(x);
    auto vars = quantifier_variable_stack.top();
    vars.insert(vars.end(), x.variables().begin(), x.variables().end());
    quantifier_variable_stack.push(vars);
    state_formula body = (*this).apply(x.body());
    quantifier_variable_stack.pop();
    state_formula result = forall(x.variables(), body);
    (*this).leave(x);
    return result;
  }

  state_formula apply(const exists& x)
  {
    (*this).enter(x);
    auto vars = quantifier_variable_stack.top();
    vars.insert(vars.end(), x.variables().begin(), x.variables().end());
    quantifier_variable_stack.push(vars);
    state_formula body = (*this).apply(x.body());
    quantifier_variable_stack.pop();
    state_formula result = exists(x.variables(), body);
    (*this).leave(x);
    return result;
  }

};

/// \endcond

/// \brief
/// \param x an object containing state formulas
state_formula transform(const modal_equation_system& x)
{
  modal_equation_transform_builder f(x);
  state_formula result = f.apply(x.initial_state());
  return result;
}

/// \brief
/// \param x an object containing state formulas
modal_equation_system transform(const state_formula& x)
{
  state_formula_transform_builder f;
  state_formula initial_state = f.apply(x);
  modal_equation_system result;
  if (is_variable(initial_state))
  {
    variable initial_variable = atermpp::down_cast<const state_formulas::variable>(initial_state);
    // equations are added in reverse order, which is repaired here:
    std::vector<modal_equation> equations = f.equations();
    std::vector<modal_equation> requations(equations.rbegin(), equations.rend());
    result = modal_equation_system(initial_variable, requations);
  }
  else
  {
    throw std::runtime_error("Top operator of the formula should be 'nu' or 'mu'.");
  }
  return result;
}


} // namespace state_formulas

} // namespace mcrl2


#endif /* MCRL2_MODAL_TRANSFORM_MODAL_EQUATION_SYSTEM_H_ */
