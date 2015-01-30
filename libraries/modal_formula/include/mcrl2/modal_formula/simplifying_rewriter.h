// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/simplifying_rewriter.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_SIMPLIFYING_REWRITER_H_
#define MCRL2_MODAL_FORMULAS_SIMPLIFYING_REWRITER_H_

#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include "boost/static_assert.hpp"
#include "mcrl2/data/bool.h"
#include "mcrl2/data/rewriter.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/modal_formula/state_formulas_optimized.h"
#include "mcrl2/utilities/detail/join.h"


namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

template <typename T>
void simplify(T& x, data::data_specification& data, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

template <typename T>
T simplify(const T& x, data::data_specification& data, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

// \brief Visitor for simplifying a state formula.
struct simplifying_rewriter_builder: public state_formula_builder<simplifying_rewriter_builder>
{
  typedef state_formula_builder<simplifying_rewriter_builder> super;
  typedef core::term_traits_optimized<state_formulas::state_formula> optimized;
  typedef core::term_traits<state_formulas::state_formula> tr;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  data::data_specification& data;
  std::map<data::data_expression, data::data_expression> rewrite_cache;

  /// \brief Rewriter for simplifying expressions
  data::rewriter::strategy rewrite_strategy;
  data::rewriter datar;

  inline state_formulas::state_formula data_to_state(const data::data_expression& x)
  {
    state_formulas::state_formula result;
    if (x == data::sort_bool::true_())
    {
      result = state_formulas::true_();
    }
    else if (x == data::sort_bool::false_())
    {
      result = state_formulas::false_();
    }
    else
    {
      result = state_formulas::state_formula(x);
    }
    return result;
  }

  size_t n_cache_calls;
  size_t n_cache_hits;

  inline data::data_expression rewrite(const data::data_expression& d)
  {
    if (d == data::sort_bool::true_() || d == data::sort_bool::false_())
    {
      return d;
    }
    n_cache_calls++;
    auto it = rewrite_cache.find(d);
    if (it == rewrite_cache.end())
    {
      data::data_expression e = datar(d);
      rewrite_cache[d] = e;
      return e;
    }
    n_cache_hits++;
    return it->second;
  }

  std::stack<std::set<core::identifier_string>> variable_replacement_candidates;
  std::stack<std::set<data::variable>> variable_bindings;
  std::map<core::identifier_string, data::assignment_list> variable_assignments;
  std::map<core::identifier_string, atermpp::function_symbol> variable_function_symbols;
  size_t variable_occurence_count;
  size_t variable_replace_count;

public:

  simplifying_rewriter_builder(data::data_specification& data_)
    : data(data_),
      rewrite_strategy(data::parse_rewrite_strategy("jitty")),
      datar(data::rewriter(data, rewrite_strategy)),
      n_cache_calls(0),
      n_cache_hits(0),
      variable_occurence_count(0),
      variable_replace_count(0)
  {
    variable_replacement_candidates.push(std::set<core::identifier_string>());
    variable_bindings.push(std::set<data::variable>());
  }

  std::string report_cache_hits()
  {
    std::stringstream s;
    s <<  "simplifying_rewriter cache: " << n_cache_hits << " / " << n_cache_calls << " hits (" <<
         ((n_cache_calls==0) ? 100 : (((float)n_cache_hits/(float)n_cache_calls)*100)) << " %)";
    return s.str();
  }

  std::string report_variable_replacements()
  {
    std::stringstream s;
    s << "variables replaced: " << variable_replace_count << " / " << variable_occurence_count;
    return s.str();
  }

  state_formulas::state_formula apply(const data::data_expression& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "data expression " << pp(x) << ":" << std::endl;
    state_formulas::state_formula result = data_to_state(rewrite(x));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::must& x)
  {
    (*this).enter(x);
    //mCRL2log(log::verbose) << "must [" << pp(x.formula()) << "] " << pp(x.operand()) << ":" << std::endl;
    variable_replacement_candidates.push(std::set<core::identifier_string>());
    state_formulas::state_formula result = optimized::must(x.formula(), (*this).apply(x.operand()));
    variable_replacement_candidates.pop();
    (*this).leave(x);
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::may& x)
  {
    (*this).enter(x);
    //mCRL2log(log::verbose) << "may <" << pp(x.formula()) << "> " << pp(x.operand()) << ":" << std::endl;
    variable_replacement_candidates.push(std::set<core::identifier_string>());
    state_formulas::state_formula result = optimized::may(x.formula(), (*this).apply(x.operand()));
    variable_replacement_candidates.pop();
    (*this).leave(x);
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::variable& x)
  {
    //mCRL2log(log::verbose) << "variable: " << pp(x) << std::endl;
    (*this).enter(x);
    std::vector<data::data_expression> arguments;
    for(auto a_it = x.arguments().begin(); a_it != x.arguments().end(); ++a_it)
    {
      arguments.push_back(rewrite(*a_it));
    }
    (*this).leave(x);
    state_formulas::state_formula result = state_formulas::variable(x.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;

    // Trying if the variable can be replaced by true or false.
    // In a formula 'nu X(d: D). phi', occurrences of 'X(d)' in phi that are not in the scope of a modal operator
    // can be replaced by 'true' (for a formula 'mu X(d: D). phi', the substitute is 'false').
    // Note that the data variable d should not be bound in phi by a quantifier 'forall d: D' or 'exists d: D'.
    if (variable_replacement_candidates.top().find(x.name())!=variable_replacement_candidates.top().end())
    {
      variable_occurence_count++;
      //mCRL2log(log::verbose) << "checking if variable " << pp(x) << " can be replaced by true or false."<< std::endl;
      bool is_self_loop = true;
      data::assignment_list assignments = variable_assignments[x.name()];
      auto it1 = assignments.begin();
      auto it2 = x.arguments().begin();
      for(; it1 != assignments.end() && it2 != x.arguments().end(); ++it1, ++it2)
      {
        data::variable lhs = it1->lhs();
        data::data_expression rhs = (*it2);
        if (data::is_variable(rhs))
        {
          data::variable var(rhs);
          if (lhs != var || variable_bindings.top().find(var) != variable_bindings.top().end())
          {
            is_self_loop = false;
            break;
          }
        }
        else
        {
          is_self_loop = false;
          break;
        }
      }
      if (is_self_loop)
      {
        variable_replace_count++;
        // applying the substitution
        result = (variable_function_symbols[x.name()] == core::detail::function_symbols::StateNu) ?
           (state_formulas::state_formula)state_formulas::true_() : (state_formulas::state_formula)state_formulas::false_();
        mCRL2log(log::verbose) << "* replacing variable '" << pp(x) << "' with '" <<
            ((variable_function_symbols[x.name()] == core::detail::function_symbols::StateNu) ?
                "true":"false") << "'" << std::endl;
      }
    }
    return result;
  }

  void enter(const state_formulas::nu& x)
  {
    variable_assignments[x.name()] = x.assignments();
    variable_function_symbols[x.name()] = core::detail::function_symbols::StateNu;
    // push a singleton set containing x
    std::set<core::identifier_string> replacement_candidates;
    replacement_candidates.insert(x.name());
    variable_replacement_candidates.push(replacement_candidates);
    auto bindings = variable_bindings.top();
    std::set<data::variable> new_bindings;
    std::set_difference(bindings.begin(), bindings.end(),
        x.assignments().begin(), x.assignments().end(),
        std::inserter(new_bindings, new_bindings.end()));
    variable_bindings.push(new_bindings);
  }

  state_formulas::state_formula apply(const state_formulas::nu& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "nu " << x.name() << ": " << pp(x) << std::endl;
    state_formulas::state_formula operand = (*this).apply(x.operand());
    std::set<core::identifier_string> occ = state_formulas::find_state_variable_names(operand);
    state_formulas::state_formula result;
    if (occ.find(x.name()) != occ.end())
    {
      // possibly rewrite 'nu Z. [a]Z' to 'nu Z. true'
      if (state_formulas::is_must(operand))
      {
        state_formulas::must m = (state_formulas::must)operand;
        if (state_formulas::is_variable(m.operand()))
        {
          state_formulas::variable var = (state_formulas::variable)m.operand();
          bool is_self_loop = true;
          data::assignment_list assignments = x.assignments();
          auto it1 = assignments.begin();
          auto it2 = var.arguments().begin();
          for(; it1 != assignments.end() && it2 != var.arguments().end(); ++it1, ++it2)
          {
            data::variable lhs = it1->lhs();
            data::data_expression rhs = (*it2);
            if (data::is_variable(rhs))
            {
              data::variable var(rhs);
              if (lhs != var)
              {
                is_self_loop = false;
                break;
              }
            }
            else
            {
              is_self_loop = false;
              break;
            }
          }
          if (is_self_loop)
          {
            // rewrite 'nu Z. [a]Z' to 'nu Z. true'
            operand = state_formulas::true_();
          }
        }
      }
      result = state_formulas::nu(x.name(), x.assignments(), operand);
    }
    else
    {
      data::mutable_map_substitution< std::map< data::variable, data::data_expression > > substitutions;
      for(const auto& a: x.assignments())
      {
        substitutions[a.lhs()] = a.rhs();
      }
      result = (*this).apply(state_formulas::replace_variables(operand, substitutions));
    }
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  void leave(const state_formulas::nu& x)
  {
    variable_bindings.pop();
    variable_replacement_candidates.pop();
  }

  void enter(const state_formulas::mu& x)
  {
    variable_assignments[x.name()] = x.assignments();
    variable_function_symbols[x.name()] = core::detail::function_symbols::StateMu;
    // push a singleton set containing x
    std::set<core::identifier_string> replacement_candidates;
    replacement_candidates.insert(x.name());
    variable_replacement_candidates.push(replacement_candidates);
    auto bindings = variable_bindings.top();
    std::set<data::variable> new_bindings;
    std::set_difference(bindings.begin(), bindings.end(),
        x.assignments().begin(), x.assignments().end(),
        std::inserter(new_bindings, new_bindings.end()));
    variable_bindings.push(new_bindings);
  }

  state_formulas::state_formula apply(const state_formulas::mu& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "mu " << x.name() << ": " << pp(x) << std::endl;
    state_formulas::state_formula operand = (*this).apply(x.operand());
    std::set<core::identifier_string> occ = state_formulas::find_state_variable_names(operand);
    state_formulas::state_formula result;
    if (occ.find(x.name()) != occ.end())
    {
      // possibly rewrite 'mu Z. <a>Z' to 'mu Z. false'
      if (state_formulas::is_may(operand))
      {
        state_formulas::must m = (state_formulas::must)operand;
        if (state_formulas::is_variable(m.operand()))
        {
          state_formulas::variable var = (state_formulas::variable)m.operand();
          bool is_self_loop = true;
          data::assignment_list assignments = x.assignments();
          auto it1 = assignments.begin();
          auto it2 = var.arguments().begin();
          for(; it1 != assignments.end() && it2 != var.arguments().end(); ++it1, ++it2)
          {
            data::variable lhs = it1->lhs();
            data::data_expression rhs = (*it2);
            if (data::is_variable(rhs))
            {
              data::variable var(rhs);
              if (lhs != var)
              {
                is_self_loop = false;
                break;
              }
            }
            else
            {
              is_self_loop = false;
              break;
            }
          }
          if (is_self_loop)
          {
            // rewrite 'mu Z. <a>Z' to 'mu Z. false'
            operand = state_formulas::false_();
          }
        }
      }
      result = state_formulas::mu(x.name(), x.assignments(), operand);
    }
    else
    {
      data::mutable_map_substitution< std::map< data::variable, data::data_expression > > substitutions;
      for(const auto& a: x.assignments())
      {
        substitutions[a.lhs()] = a.rhs();
      }
      result = (*this).apply(state_formulas::replace_variables(operand, substitutions));
    }
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  void leave(const state_formulas::mu& x)
  {
    variable_bindings.pop();
    variable_replacement_candidates.pop();
  }

  state_formulas::state_formula apply(const state_formulas::not_& x)
  {
    (*this).enter(x);
    state_formulas::state_formula result = optimized::not_((*this).apply(x.operand()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::and_& x)
  {
    (*this).enter(x);
    //mCRL2log(log::verbose) << "and: " << pp(x) << std::endl;
    state_formulas::state_formula result = optimized::and_((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::or_& x)
  {
    //mCRL2log(log::verbose) << "or: " << pp(x) << std::endl;
    (*this).enter(x);
    state_formulas::state_formula result = optimized::or_((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::imp& x)
  {
    //mCRL2log(log::verbose) << "imp: " << pp(x) << std::endl;
    (*this).enter(x);
    state_formulas::state_formula result = optimized::imp((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    return result;
  }

  void enter(const state_formulas::forall& x)
  {
    auto bindings = variable_bindings.top();
    bindings.insert(x.variables().begin(), x.variables().end());
    variable_bindings.push(bindings);
  }

  state_formulas::state_formula apply(const state_formulas::forall& x)
  {
    //mCRL2log(log::verbose) << "forall: " << pp(x) << std::endl;
    (*this).enter(x);
    state_formulas::state_formula result = optimized::forall(x.variables(), (*this).apply(x.body()));
    (*this).leave(x);
    return result;
  }

  void leave(const state_formulas::forall& x)
  {
    variable_bindings.pop();
  }

  void enter(const state_formulas::exists& x)
  {
    auto bindings = variable_bindings.top();
    bindings.insert(x.variables().begin(), x.variables().end());
    variable_bindings.push(bindings);
  }

  state_formulas::state_formula apply(const state_formulas::exists& x)
  {
    //mCRL2log(log::verbose) << "exists: " << pp(x) << std::endl;
    (*this).enter(x);
    state_formulas::state_formula result = optimized::exists(x.variables(), (*this).apply(x.body()));
    (*this).leave(x);
    return result;
  }

  void leave(const state_formulas::exists& x)
  {
    variable_bindings.pop();
  }

  void apply(modal_equation& x)
  {
    state_formula formula = (*this).apply(x.formula());
    x.set_formula(formula);
  }

  void apply(modal_equation_system& x)
  {
    std::vector<modal_equation> equations;
    for(auto e: x.equations())
    {
      (*this).apply(e);
      equations.push_back(e);
    }
    x.set_equations(equations);
  }

};
/// \endcond

/// \brief
/// \param x an object containing state formulas
template <typename T>
void simplify(T& x, data::data_specification& data, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  simplifying_rewriter_builder f(data);
  f.apply(x);
}

/// \brief
/// \param x an object containing state formulas
template <typename T>
T simplify(const T& x, data::data_specification& data, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  simplifying_rewriter_builder f(data);
  T result = f.apply(x);
  mCRL2log(log::verbose) << "(" << f.report_cache_hits() << ")" << std::endl;
  mCRL2log(log::verbose) << "(" << f.report_variable_replacements() << ")" << std::endl;
  return result;
}


} // namespace state_formulas

} // namespace mcrl2

#endif /* MCRL2_MODAL_FORMULAS_SIMPLIFYING_REWRITER_H_ */
