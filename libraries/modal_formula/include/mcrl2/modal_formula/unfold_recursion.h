// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/unfold_recursion.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_UNFOLD_RECURSION_H_
#define MCRL2_MODAL_FORMULAS_UNFOLD_RECURSION_H_

#include "boost/static_assert.hpp"
#include "mcrl2/data/bool.h"
#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/formula_size.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/rewrite_quantifiers.h"
#include "mcrl2/modal_formula/simplifying_rewriter.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/modal_formula/state_formulas_optimized.h"
#include "mcrl2/utilities/detail/join.h"


namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

// \brief Visitor for finding equations in a state formula.
struct state_formula_equation_traverser: public state_formula_traverser<state_formula_equation_traverser>
{
  typedef state_formula_traverser<state_formula_equation_traverser> super;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  std::map<core::identifier_string, std::pair<data::variable_list, state_formulas::state_formula > > equations;

public:
  state_formula_equation_traverser()
  {
  }

  std::map<core::identifier_string, std::pair<data::variable_list, state_formulas::state_formula > >& get_equations()
  {
    return equations;
  }

  void apply(const state_formulas::nu& x)
  {
    (*this).enter(x);
    std::vector<data::variable> variables;
    for (const auto& a: x.assignments()) {
      variables.push_back(a.lhs());
    }
    equations[x.name()] = std::make_pair(data::variable_list(variables.begin(), variables.end()), x.operand());
    (*this).leave(x);
  }

  void apply(const state_formulas::mu& x)
  {
    (*this).enter(x);
    std::vector<data::variable> variables;
    for (const auto& a: x.assignments()) {
      variables.push_back(a.lhs());
    }
    equations[x.name()] = std::make_pair(data::variable_list(variables.begin(), variables.end()), x.operand());
    (*this).leave(x);
  }

};

template <typename T>
void unfold_recursion(T& x, data::data_specification& data,
    typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

template <typename T>
T unfold_recursion(const T& x, data::data_specification& data,
    typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

// \brief Visitor for unfolding unguarded recursion in a state formula.
struct unfold_recursion_builder: public state_formula_builder<unfold_recursion_builder>
{
  typedef state_formula_builder<unfold_recursion_builder> super;
  typedef core::term_traits_optimized<state_formulas::state_formula> optimized;
  typedef core::term_traits<state_formulas::state_formula> tr;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  /// Data specification
  data::data_specification& m_data;
  // map from variable name to corresponding parameters and right hand side formula
  std::map<core::identifier_string, std::pair<data::variable_list, state_formulas::state_formula > >& m_equations;
  // Rewrite cache
  std::map<state_formulas::state_formula,
    std::map< std::map<data::variable, data::data_expression>,
      state_formulas::state_formula > > rewrite_cache;
  bool use_rewrite_cache;
  size_t counter;
  simplifying_rewriter_builder simplify;
  rewrite_quantifiers_builder rewrite_quantifiers;

  state_formulas::state_formula replace_variables(const state_formulas::state_formula& x, const std::map<data::variable, data::data_expression>& s)
  {
    data::mutable_map_substitution< std::map< data::variable, data::data_expression > > sigma(s);
    if (tr::is_true(x) || tr::is_false(x))
    {
      return x;
    }
    else
    {
      if (use_rewrite_cache)
      {
        auto it1 = rewrite_cache.find(x);
        if (it1 == rewrite_cache.end())
        {
          std::map< std::map< data::variable, data::data_expression >,
                state_formulas::state_formula > entry;
          state_formulas::state_formula r = state_formulas::replace_variables(x, sigma);
          entry[s] = r;
          rewrite_cache[x] = entry;
          return r;
        }
        else
        {
          auto entry = it1->second;
          auto it2 = entry.find(s);
          if (it2 == entry.end())
          {
            state_formulas::state_formula r = state_formulas::replace_variables(x, sigma);
            entry[s] = r;
            return r;
          }
          else
          {
            return it2->second;
          }
        }
      }
      else
      {
        return state_formulas::replace_variables(x, sigma);
      }
    }
  }

public:
  unfold_recursion_builder(
      data::data_specification& data,
      std::map<core::identifier_string, std::pair<data::variable_list, state_formulas::state_formula > >& equations) :
    m_data(data),
    m_equations(equations),
    use_rewrite_cache(true),
    counter(0),
    simplify(m_data),
    rewrite_quantifiers(m_data)
  {}

  size_t get_counter()
  {
    return counter;
  }

  state_formulas::state_formula apply(const state_formulas::must& x)
  {
    (*this).enter(x);
    // skip
    (*this).leave(x);
    return x;
  }

  state_formulas::state_formula apply(const state_formulas::may& x)
  {
    (*this).enter(x);
    // skip
    (*this).leave(x);
    return x;
  }

  state_formulas::state_formula apply(const state_formulas::variable& x)
  {
    (*this).enter(x);
    //mCRL2log(log::verbose) << "replacing variable with its formula..." << std::endl;
    auto entry = m_equations.find(x.name());
    if (entry == m_equations.end())
    {
      throw std::runtime_error("Cannot find equation for variable " + std::string(x.name()));
    }
    auto eq = entry->second;
    data::variable_list a = eq.first;
    auto it1 = a.begin();
    auto it2 = x.arguments().begin();
    std::map< data::variable, data::data_expression > sigma;
    for(; it1 != a.end() && it2 != x.arguments().end(); ++it1, ++it2)
    {
      if ((*it1) != (*it2))
      {
        sigma[*it1] = (*it2);
      }
    }
    state_formulas::state_formula formula;
    //mCRL2log(log::verbose) << "applying substitution." << std::endl;
    formula = replace_variables(eq.second, sigma);
    //size_t size1 = formula_size(formula);
    //mCRL2log(log::verbose) << "size of resulting subformula: " << size1 << std::endl;
    //mCRL2log(log::verbose) << "rewriting formula..." << std::endl;
    formula = rewrite_quantifiers.apply(formula);
    //mCRL2log(log::verbose) << "simplifying formula..." << std::endl;
    formula = simplify.apply(formula);
    //size_t size2 = formula_size(formula);
    //mCRL2log(log::verbose) << "size after simplification: " << size2
    //    << " (reduced to: " << ((size1==0) ? 0 : (((float)size2/(float)size1)*100)) << " %)"<< std::endl;
    counter++;
    (*this).leave(x);
    //if (size2 < 10 || size1/size2 > 5)
    //{
    //  mCRL2log(log::verbose) << "unfold variable..." <<
    //      " (formula reduced to " << ((size1==0) ? 0 : (((float)size2/(float)size1)*100)) << " %)"<< std::endl;
    //}
    return formula;
  }

};
/// \endcond

/// \brief
/// \param x an object containing state formulas
template <typename T>
void unfold_recursion(T& x, data::data_specification& data, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  size_t size1 = formula_size(x);
  state_formula_equation_traverser eq;
  eq.apply(x);
  unfold_recursion_builder f(data, eq.get_equations());
  x = f.apply(x);
  mCRL2log(log::verbose) << "unfolded " << f.get_counter() << " variables." << std::endl;
  mCRL2log(log::verbose) << "(" << f.simplify.report_cache_hits() << ")" << std::endl;
  mCRL2log(log::verbose) << "(" << f.rewrite_quantifiers.report_cache_hits() << ")" << std::endl;
  size_t size2 = formula_size(x);
  mCRL2log(log::verbose) << "(formula increased to " << ((size1==0) ? 0 : (((float)size2/(float)size1)*100)) << " %)"<< std::endl;
}

/// \brief
/// \param x an object containing state formulas
void unfold_recursion(modal_equation_system& x, data::data_specification& data)
{
  size_t size1 = formula_size(x);
  std::map<core::identifier_string, std::pair<data::variable_list, state_formulas::state_formula > > equations_map;
  for (const auto& entry: x.equation_map())
  {
    equations_map[entry.first] = std::make_pair(entry.second.parameters(), entry.second.formula());
  }
  unfold_recursion_builder f(data, equations_map);
  auto equations = x.equations();
  for (auto& eq: equations)
  {
    eq.formula() = f.apply(eq.formula());
  }
  x.set_equations(equations);
  mCRL2log(log::verbose) << "unfolded " << f.get_counter() << " variables." << std::endl;
  mCRL2log(log::verbose) << "(" << f.simplify.report_cache_hits() << ")" << std::endl;
  mCRL2log(log::verbose) << "(" << f.rewrite_quantifiers.report_cache_hits() << ")" << std::endl;
  size_t size2 = formula_size(x);
  mCRL2log(log::verbose) << "(formula increased to " << ((size1==0) ? 0 : (((float)size2/(float)size1)*100)) << " %)"<< std::endl;
}


/// \brief
/// \param x an object containing state formulas
template <typename T>
T unfold_recursion(const T& x, data::data_specification& data, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  size_t size1 = formula_size(x);
  state_formula_equation_traverser eq;
  eq.apply(x);
  unfold_recursion_builder f(data, eq.get_equations());
  T result = f.apply(x);
  mCRL2log(log::verbose) << "unfolded " << f.get_counter() << " variables." << std::endl;
  mCRL2log(log::verbose) << "(" << f.simplify.report_cache_hits() << ")" << std::endl;
  mCRL2log(log::verbose) << "(" << f.rewrite_quantifiers.report_cache_hits() << ")" << std::endl;
  size_t size2 = formula_size(result);
  mCRL2log(log::verbose) << "(formula increased to " << ((size1==0) ? 0 : (((float)size2/(float)size1)*100)) << " %)"<< std::endl;
  return result;
}


} // namespace state_formulas

} // namespace mcrl2


#endif /* MCRL2_MODAL_FORMULAS_UNFOLD_RECURSION_H_ */
