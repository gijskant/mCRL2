// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/rewrite_quantifiers.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_REWRITE_QUANTIFIERS_H_
#define MCRL2_MODAL_FORMULAS_REWRITE_QUANTIFIERS_H_

#include "boost/static_assert.hpp"
#include "mcrl2/data/bool.h"
#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/simplifying_rewriter.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/normalize.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/modal_formula/state_formulas_optimized.h"
#include "mcrl2/utilities/detail/join.h"


namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

struct rewrite_quantifiers_accessors
{
  typedef core::term_traits<state_formulas::state_formula> tr;

  static inline
  bool is_and(const state_formulas::state_formula& x)
  {
    return tr::is_and(x)
    || (tr::is_data(x)
        && data::sort_bool::is_and_application(atermpp::down_cast<const data::application>(x)));
  }

  static inline
  bool is_or(const state_formulas::state_formula& x)
  {
    return tr::is_or(x) || tr::is_imp(x)
    || (tr::is_data(x)
        && (data::sort_bool::is_or_application(atermpp::down_cast<const data::application>(x))
            || data::sort_bool::is_implies_application(atermpp::down_cast<const data::application>(x))));
  }

  static inline
  state_formulas::state_formula left(const state_formulas::state_formula& x)
  {
    assert(is_and(x) || is_or(x));
    if (tr::is_or(x) || tr::is_and(x))
    {
      return tr::left(x);
    }
    else if (tr::is_imp(x))
    {
      return tr::not_(tr::left(x));
    }
    else {
      assert(tr::is_data(x));
      const data::application& appl = atermpp::down_cast<const data::application>(x);
      if (data::sort_bool::is_implies_application(appl))
      {
        return tr::not_(appl[0]);
      }
      else
      {
        return appl[0];
      }
    }
  }

  static inline
  state_formulas::state_formula right(const state_formulas::state_formula& x)
  {
    assert(is_and(x) || is_or(x));
    if (tr::is_or(x) || tr::is_and(x) || tr::is_imp(x))
    {
      return tr::right(x);
    }
    else {
      assert(tr::is_data(x));
      const data::application& appl = atermpp::down_cast<const data::application>(x);
      return appl[1];
    }
  }

  static inline
  std::vector<state_formulas::state_formula> split_conjuncts(const state_formulas::state_formula& x)
  {
    std::vector<state_formulas::state_formula> result;
    utilities::detail::split(x, std::back_insert_iterator<std::vector<state_formulas::state_formula> >(result), is_and, left, right);
    return result;
  }

  static inline
  std::vector<state_formulas::state_formula> split_disjuncts(const state_formulas::state_formula& x)
  {
    std::vector<state_formulas::state_formula> result;
    utilities::detail::split(x, std::back_insert_iterator<std::vector<state_formulas::state_formula> >(result), is_or, left, right);
    return result;
  }

  static inline std::pair<bool, data::data_expression> find_equality_for_variable(
      data::data_expression expr1, data::data_expression expr2, data::variable d)
  {
    //mCRL2log(log::debug) << "find_equality_for_variable for variable " << pp(d) << ": expr1 = " << pp(expr1)
    //    << ", expr2 = " << pp(expr2) << std::endl;
    bool result = false;
    data::data_expression expr = data::sort_bool::false_();
    if (data::is_variable(expr1))
    {
      data::variable var = atermpp::down_cast<const data::variable>(expr1);
      if (var==d)
      {
        result = true;
        expr = expr2;
      }
      else if (data::is_variable(expr2))
      {
        var = atermpp::down_cast<const data::variable>(expr2);
        if (var==d)
        {
          result = true;
          expr = expr1;
        }
      }
    }
    //mCRL2log(log::debug) << "find_equality_for_variable: result = " << result << ", expr = " << pp(expr) << std::endl;
    return std::make_pair(result, expr);
  }

  static std::pair<bool, data::data_expression> find_point(const state_formulas::state_formula& x, data::variable d, bool negative = false)
  {
    //mCRL2log(log::verbose) << " - find_point for variable " << pp(d) << ": " << pp(x) << std::endl;
    //mCRL2log(log::verbose) << "find_point for variable: negative = " << negative << std::endl;
    bool is_point = false;
    data::data_expression expr = data::sort_bool::false_();
    if (tr::is_not(x))
    {
      return find_point(tr::not_arg(x), d, !negative);
    }
    else if (tr::is_data(x))
    {
      //mCRL2log(log::verbose) << "Considering data expression: " << pp(x) << std::endl;
      const data::data_expression& e = atermpp::down_cast<const data::data_expression>(x);
      if (data::sort_bool::is_not_application(e))
      {
        return find_point(data::sort_bool::arg(e), d, !negative);
      }
      else if ((negative && data::is_not_equal_to_application(e)) || (!negative && data::is_equal_to_application(e)))
      {
        const data::application& appl = atermpp::down_cast<const data::application>(e);
        if (appl.size()==2)
        {
          data::data_expression expr1 = appl[0];
          data::data_expression expr2 = appl[1];
          auto p = find_equality_for_variable(expr1, expr2, d);
          is_point = p.first;
          if (is_point)
          {
            expr = p.second;
          }
          else if (!is_point)
          {
            p = find_equality_for_variable(expr2, expr1, d);
            is_point = p.first;
            if (is_point)
            {
              expr = p.second;
            }
          }
          if (is_point)
          {
            //mCRL2log(log::verbose) << " - Found point for variable " << pp(d) << ": " << pp(expr) << std::endl;
            std::set<data::variable> free_vars = data::find_free_variables(expr);
            if (free_vars.find(d) != free_vars.end())
            {
              //mCRL2log(log::verbose) << " - Expression contains variable " << pp(d) << ". Cannot remove quantifier." << std::endl;
              is_point = false;
            }
            else
            {
              //mCRL2log(log::verbose) << " - Will substitute variable " << pp(d) << " with expression." << std::endl;
            }
          }
        }
      }
      else if (is_and(x) && !negative)
      {
        auto p = find_point(left(x), d, negative);
        if (p.first)
        {
          return p;
        }
        else
        {
          return find_point(right(x), d, negative);
        }
      }
      else if (is_or(x) && negative)
      {
        auto p = find_point(left(x), d, negative);
        if (p.first)
        {
          return p;
        }
        else
        {
          return find_point(right(x), d, negative);
        }
      }
    }
    return std::make_pair(is_point, expr);
  }

};

template <typename T>
void discover_parameters(T& x, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

template <typename T>
T discover_parameters(const T& x, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

// \brief Visitor for discovering variable parameters in a state formula.
struct discover_parameters_builder: public state_formula_builder<discover_parameters_builder>
{
  typedef state_formula_builder<discover_parameters_builder> super;
  typedef core::term_traits_optimized<state_formulas::state_formula> optimized;
  typedef core::term_traits<state_formulas::state_formula> tr;
  typedef rewrite_quantifiers_accessors accessors;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  std::stack<std::set<data::variable>> parameters_in_scope;
  size_t n_discoveries;

public:
  discover_parameters_builder() :
    n_discoveries(0)
  {
    parameters_in_scope.push(std::set<data::variable>());
  }

  std::string report()
  {
    std::stringstream s;
    s <<  "discover_parameters: " << n_discoveries << " replacements.";
    return s.str();
  }

  // Rewriting formulae of the shape 'd==5 => phi' to 'd==5 => phi[5:=d]' if the formula
  // occurs as a subformula of a formula 'nu X(d: D) . psi'. Only if the parameter 'd' is not
  // overwritten/bound by a quantifier.
  state_formulas::state_formula apply(const state_formulas::or_& x)
  {
    //mCRL2log(log::verbose) << "or_: " << pp(x) << std::endl;
    (*this).enter(x);
    data::mutable_map_substitution< std::map< data::data_expression, data::data_expression > > substitutions;
    bool apply = false;
    bool changed = false;
    state_formulas::state_formula normalized_x = normalize((state_formulas::state_formula)x);
    std::vector<state_formulas::state_formula> parts = accessors::split_disjuncts(normalized_x);
    std::vector<state_formulas::state_formula> result_parts;
    std::vector<state_formulas::state_formula> point_parts;
    for(const auto& p: parts)
    {
      //mCRL2log(log::verbose) << " *  Looking for an equality in: " << pp(p) << std::endl;
      bool is_point = false;
      for(const auto& v: parameters_in_scope.top())
      {
        //mCRL2log(log::verbose) << "- Considering parameter: " << pp(v) << std::endl;
        data::data_expression expr = data::sort_bool::false_();
        auto point = accessors::find_point(p, v, true);
        is_point = point.first;
        if (is_point)
        {
          expr = point.second;
          //mCRL2log(log::verbose) << "Substitute expression " << pp(expr) << " with parameter: " << pp(v) << std::endl;
          substitutions[expr] = v;
          apply = true;
          break;
        }
      }
      if (is_point)
      {
        point_parts.push_back(p);
      }
      else
      {
        state_formulas::state_formula p_result = (*this).apply(p);
        if (p_result != p)
        {
          changed = true;
        }
        result_parts.push_back(p_result);
      }
    }
    state_formulas::state_formula result;
    // substitution [e/v]
    if (apply || changed)
    {
      //mCRL2log(log::verbose) << "Applying substitution..." << std::endl;
      result = optimized::join_or(result_parts.begin(), result_parts.end());
      if (apply)
      {
        result = state_formulas::replace_data_expressions(result, substitutions, false);
        result = optimized::or_(optimized::join_or(point_parts.begin(), point_parts.end()), result);
      }
      if (result != x)
      {
        n_discoveries++;
        //if (mCRL2logEnabled(log::verbose))
        //{
        //  mCRL2log(log::verbose) << "  * expr: " << pp(x) << std::endl;
        //  mCRL2log(log::verbose) << "  * ----> " << pp(result) << std::endl;
        //}
      }
    }
    else
    {
      //mCRL2log(log::verbose) << "Not applying substitution." << std::endl;
      result = x;
    }
    (*this).leave(x);
    return result;
  }

  // Rewriting formulae of the shape 'd==5 && phi' to 'd==5 && phi[5:=d]' if the formula
  // occurs as a subformula of a formula 'nu X(d: D) . psi'. Only if the parameter 'd' is not
  // overwritten/bound by a quantifier.
  state_formulas::state_formula operator()(const state_formulas::and_& x)
  {
    //mCRL2log(log::verbose) << "and_: " << pp(x) << std::endl;
    (*this).enter(x);
    data::mutable_map_substitution< std::map< data::data_expression, data::data_expression > > substitutions;
    bool apply = false;
    bool changed = false;
    state_formulas::state_formula normalized_x = normalize((state_formulas::state_formula)x);
    std::vector<state_formulas::state_formula> parts = accessors::split_conjuncts(normalized_x);
    std::vector<state_formulas::state_formula> result_parts;
    std::vector<state_formulas::state_formula> point_parts;
    for(const auto& p: parts)
    {
      //mCRL2log(log::verbose) << " *  Looking for an equality in: " << pp(p) << std::endl;
      bool is_point = false;
      for(const auto& v: parameters_in_scope.top())
      {
        //mCRL2log(log::verbose) << "- Considering parameter: " << pp(v) << std::endl;
        data::data_expression expr = data::sort_bool::false_();
        auto point = accessors::find_point(p, v, false);
        is_point = point.first;
        if (is_point)
        {
          expr = point.second;
          //mCRL2log(log::verbose) << "Substitute expression " << pp(expr) << " with parameter: " << pp(v) << std::endl;
          substitutions[expr] = v;
          apply = true;
          break;
        }
      }
      if (is_point)
      {
        point_parts.push_back(p);
      }
      else
      {
        state_formulas::state_formula p_result = (*this).apply(p);
        if (p_result != p)
        {
          changed = true;
        }
        result_parts.push_back(p_result);
      }
    }
    state_formulas::state_formula result;
    // substitution [e/v]
    if (apply || changed)
    {
      //mCRL2log(log::verbose) << "Applying substitution..." << std::endl;
      result = optimized::join_and(result_parts.begin(), result_parts.end());
      if (apply)
      {
        result = state_formulas::replace_data_expressions(result, substitutions, false);
        result = optimized::and_(optimized::join_and(point_parts.begin(), point_parts.end()), result);
      }
      if (result != x)
      {
        n_discoveries++;
        //if (mCRL2logEnabled(log::verbose))
        //{
        //  mCRL2log(log::verbose) << "  * expr: " << pp(x) << std::endl;
        //  mCRL2log(log::verbose) << "  * ----> " << pp(result) << std::endl;
        //}
      }
    }
    else
    {
      //mCRL2log(log::verbose) << "Not applying substitution." << std::endl;
      result = x;
    }
    (*this).leave(x);
    return result;
  }

  void enter(const state_formulas::forall& x)
  {
    auto bindings = parameters_in_scope.top();
    std::set<data::variable> new_bindings;
    std::set_difference(bindings.begin(), bindings.end(),
        x.variables().begin(), x.variables().end(),
        std::inserter(new_bindings, new_bindings.end()));
    parameters_in_scope.push(new_bindings);
  }

  void leave(const state_formulas::forall& x)
  {
    parameters_in_scope.pop();
  }

  void enter(const state_formulas::exists& x)
  {
    auto bindings = parameters_in_scope.top();
    std::set<data::variable> new_bindings;
    std::set_difference(bindings.begin(), bindings.end(),
        x.variables().begin(), x.variables().end(),
        std::inserter(new_bindings, new_bindings.end()));
    parameters_in_scope.push(new_bindings);
  }

  void leave(const state_formulas::exists& x)
  {
    parameters_in_scope.pop();
  }

  void enter(const state_formulas::nu& x)
  {
    auto bindings = parameters_in_scope.top();
    for(const auto& a: x.assignments())
    {
      bindings.insert(a.lhs());
    }
    parameters_in_scope.push(bindings);
  }

  void leave(const state_formulas::nu& x)
  {
    parameters_in_scope.pop();
  }

  void enter(const state_formulas::mu& x)
  {
    auto bindings = parameters_in_scope.top();
    for(const auto& a: x.assignments())
    {
      bindings.insert(a.lhs());
    }
    parameters_in_scope.push(bindings);
  }

  void leave(const state_formulas::mu& x)
  {
    parameters_in_scope.pop();
  }


};
/// \endcond

/// \brief
/// \param x an object containing state formulas
template <typename T>
void discover_parameters(T& x, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  discover_parameters_builder f;
  f(x);
}

/// \brief
/// \param x an object containing state formulas
template <typename T>
T discover_parameters(const T& x, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  discover_parameters_builder f;
  T result = f(x);
  mCRL2log(log::verbose) << "(" << f.report() << ")" << std::endl;
  return result;
}


/// \cond INTERNAL_DOCS

template <typename T>
void rewrite_quantifiers(T& x, data::data_specification& data, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

template <typename T>
T rewrite_quantifiers(const T& x, data::data_specification& data, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

// \brief Visitor for rewriting quantifiers in a state formula.
struct rewrite_quantifiers_builder: public state_formula_builder<rewrite_quantifiers_builder>
{
  typedef state_formula_builder<rewrite_quantifiers_builder> super;
  typedef core::term_traits_optimized<state_formulas::state_formula> optimized;
  typedef core::term_traits<state_formulas::state_formula> tr;
  typedef rewrite_quantifiers_accessors accessors;
  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  /// \brief Rewriter for simplifying expressions
  state_formulas::simplifying_rewriter_builder formular;

  std::map<state_formulas::state_formula, state_formulas::state_formula> rewrite_cache;

  size_t n_cache_calls;
  size_t n_cache_hits;

  inline state_formulas::state_formula rewrite_formula(const state_formulas::state_formula& x)
  {
    if (tr::is_true(x) || tr::is_false(x))
    {
      return x;
    }
    n_cache_calls++;
    auto it = rewrite_cache.find(x);
    if (it == rewrite_cache.end())
    {
      state_formulas::state_formula r = formular.apply(x);
      rewrite_cache[x] = r;
      return r;
    }
    n_cache_hits++;
    return it->second;
  }

public:
  rewrite_quantifiers_builder(data::data_specification& data) :
    formular(data),
    n_cache_calls(0),
    n_cache_hits(0)
  {}

  std::string report_cache_hits()
  {
    std::stringstream s;
    s <<  "rewrite_quantifiers cache: " << n_cache_hits << " / " << n_cache_calls << " hits (" <<
         ((n_cache_calls==0) ? 100 : (((float)n_cache_hits/(float)n_cache_calls)*100)) << " %)";
    return s.str();
  }

  // Rewrite formulas of the shape forall x: D . x != v || f to f[x := v] if v is of type D. (one point rule)
  // forall x: D . x != v || f
  // = true && true && ... && f[x := v] && ... && true
  // = f[x := v]
  // forall x: D . x != v || !g || f
  // forall x: D . !(x == v && g) || f
  // = true && true && ... && (!g || f[x := v]) && ... && true
  // = f[x := v]
  state_formulas::state_formula apply(const state_formulas::forall& x)
  {
    //mCRL2log(log::verbose) << std::endl << "forall: " << pp(x) << std::endl;
    (*this).enter(x);
    if (x.variables().empty())
    {
      return (*this).apply(x.body());
    }
    std::vector<data::variable> x_variables;
    x_variables.push_back(x.variables().front()); // head of the list
    std::vector<data::variable> body_variables(++(x.variables().begin()), x.variables().end()); // tail of the list
    state_formulas::state_formula body = optimized::forall(tr::variable_sequence_type(body_variables.begin(), body_variables.end()), x.body());
    body = (*this).apply(body);
    state_formulas::state_formula normalized_body = normalize(body);

    // Because 'forall d. (phi || psi)' <==> 'phi || forall d. psi' if d not in free(phi), we first split disjuncts
    state_formulas::state_formula relevant_body;
    state_formulas::state_formula free_body = state_formulas::false_();
    {
      std::vector<state_formulas::state_formula> disjuncts = accessors::split_disjuncts(normalized_body);
      std::vector<state_formulas::state_formula> relevant_disjuncts;
      std::vector<state_formulas::state_formula> free_disjuncts;
      for (const auto& disjunct: disjuncts)
      {
        std::set<data::variable> free_vars = state_formulas::find_free_variables(disjunct);
        std::set<data::variable> vars;
        std::set_intersection(free_vars.begin(), free_vars.end(),
            x_variables.begin(), x_variables.end(),
            std::inserter(vars, vars.end()));
        if (vars.empty()) // d not in free(phi)
        {
          free_disjuncts.push_back((*this).apply(disjunct));
        }
        else
        {
          relevant_disjuncts.push_back(disjunct);
        }
      }
      if (free_disjuncts.empty())
      {
        relevant_body = normalized_body;
      }
      else
      {
        relevant_body = optimized::join_or(relevant_disjuncts.begin(), relevant_disjuncts.end());
        free_body = optimized::join_or(free_disjuncts.begin(), free_disjuncts.end());
      }
    }
    //mCRL2log(log::verbose) << "- free: " << pp(free_body) << std::endl;
    //mCRL2log(log::verbose) << "- relevant: " << pp(relevant_body) << std::endl;

    // Because forall d: (phi && psi) <==> forall d. phi && forall d. psi, we then split conjuncts.
    std::vector<state_formulas::state_formula> conjuncts = accessors::split_conjuncts(relevant_body);
    //conjuncts.push_back(normalized_body);
    std::vector<state_formulas::state_formula> result_conjuncts;
    for (auto c_it = conjuncts.begin(); c_it != conjuncts.end(); ++c_it)
    {
      state_formulas::state_formula conjunct = (*this).apply(*c_it);
      std::vector<state_formulas::state_formula> parts = accessors::split_disjuncts(conjunct);
      std::set<data::variable> substitution_vars;
      bool apply = false;
      data::mutable_map_substitution< std::map< data::variable, data::data_expression > > substitutions;
      for(auto v_it = x_variables.begin(); v_it != x_variables.end(); ++v_it)
      {
        bool is_point = false;
        data::data_expression expr = data::sort_bool::false_();
        for(auto p_it = parts.begin(); p_it != parts.end(); ++p_it)
        {
          auto p = accessors::find_point((*p_it), (*v_it), /* negative = */ true);
          is_point = p.first;
          if (is_point)
          {
            expr = p.second;
            //mCRL2log(log::debug) << "Substitute " << pp(*v_it) << " with expression: " << pp(expr) << std::endl;
            substitutions[(*v_it)] = expr;
            substitution_vars.insert(*v_it);
            apply = true;
          }
        }
      }
      // substitution [e/v]
      if (apply)
      {
        conjunct = state_formulas::replace_variables(conjunct, substitutions);
      }
      //mCRL2log(log::debug) << "rewriting..." << std::endl;
      conjunct = rewrite_formula(conjunct);
      //mCRL2log(log::debug) << "rewritten conjunct: " << pp(conjunct) << std::endl;
      // remove substituted variables from the set of declared variables
      std::set<data::variable> decl_vars;
      std::set_difference(x_variables.begin(), x_variables.end(),
          substitution_vars.begin(), substitution_vars.end(),
          std::inserter(decl_vars, decl_vars.end()));

      // remove unused variables (e.g., forall x. d => 3 ---> d => 3)
      std::set<data::variable> used_vars = state_formulas::find_free_variables(conjunct);
      std::set<data::variable> vars;
      std::set_intersection(decl_vars.begin(), decl_vars.end(),
          used_vars.begin(), used_vars.end(),
          std::inserter(vars, vars.end()));

      data::variable_list var_list(vars.begin(), vars.end());
      conjunct = optimized::forall(var_list, conjunct);
      if (tr::is_false(conjunct))
      {
        //if (mCRL2logEnabled(log::verbose) && x != free_body)
        //{
        //  mCRL2log(log::verbose) << "forall: " << pp(x) << std::endl;
        //  mCRL2log(log::verbose) << "  ----> " << pp(free_body) << std::endl;
        //}
        return free_body; //tr::false_();
      }
      result_conjuncts.push_back(conjunct);
      //if (mCRL2logEnabled(log::verbose) && (*c_it) != conjunct)
      //{
      //  mCRL2log(log::verbose) << "  * conjunct: " << pp(*c_it) << std::endl;
      //  mCRL2log(log::verbose) << "  *     ----> " << pp(conjunct) << std::endl;
      //}
    }
    (*this).leave(x);
    state_formulas::state_formula result =
        optimized::or_(
            free_body,
            optimized::join_and(result_conjuncts.begin(), result_conjuncts.end()));
    //if (mCRL2logEnabled(log::verbose) && x != result)
    //{
    //  mCRL2log(log::verbose) << "forall: " << pp(x) << std::endl;
    //  mCRL2log(log::verbose) << "  ----> " << pp(result) << std::endl;
    //}
    return result;
  }

  // Rewrite formulas of the shape exists x: D . x == v && f to f[x := v] if v is of type D. (one point rule)
  // exists x: D . x == v && f
  // = false || false || ... || f[x := v] || ... || false
  // = f[x := v]
  state_formulas::state_formula apply(const state_formulas::exists& x)
  {
    //mCRL2log(log::verbose) << std::endl << "exists: " << pp(x) << std::endl;
    (*this).enter(x);
    if (x.variables().empty())
    {
      return (*this).apply(x.body());
    }
    std::vector<data::variable> x_variables;
    x_variables.push_back(x.variables().front()); // head of the list
    std::vector<data::variable> body_variables(++(x.variables().begin()), x.variables().end()); // tail of the list
    state_formulas::state_formula body = optimized::exists(tr::variable_sequence_type(body_variables.begin(), body_variables.end()), x.body());
    body = (*this).apply(body);
    state_formulas::state_formula normalized_body = normalize(body);

    // Because 'exists d. (phi && psi)' <==> 'phi && exists d. psi' if d not in free(phi), we first split conjuncts
    state_formulas::state_formula relevant_body;
    state_formulas::state_formula free_body = state_formulas::true_();
    {
      std::vector<state_formulas::state_formula> conjuncts = accessors::split_conjuncts(normalized_body);
      std::vector<state_formulas::state_formula> relevant_conjuncts;
      std::vector<state_formulas::state_formula> free_conjuncts;
      for (const auto& conjunct: conjuncts)
      {
        std::set<data::variable> free_vars = state_formulas::find_free_variables(conjunct);
        std::set<data::variable> vars;
        std::set_intersection(free_vars.begin(), free_vars.end(),
            x_variables.begin(), x_variables.end(),
            std::inserter(vars, vars.end()));
        if (vars.empty()) // d not in free(phi)
        {
          free_conjuncts.push_back((*this).apply(conjunct));
        }
        else
        {
          relevant_conjuncts.push_back(conjunct);
        }
      }
      if (free_conjuncts.empty())
      {
        relevant_body = normalized_body;
      }
      else
      {
        relevant_body = optimized::join_and(relevant_conjuncts.begin(), relevant_conjuncts.end());
        free_body = optimized::join_and(free_conjuncts.begin(), free_conjuncts.end());
      }
    }
    //mCRL2log(log::verbose) << "- free: " << pp(free_body) << std::endl;
    //mCRL2log(log::verbose) << "- relevant: " << pp(relevant_body) << std::endl;

    // Because exists d. (phi || psi) <==> exists d. phi || exists d. psi, we then split disjuncts.
    std::vector<state_formulas::state_formula> disjuncts = accessors::split_disjuncts(relevant_body);
    //disjuncts.push_back(normalized_body);
    std::vector<state_formulas::state_formula> result_disjuncts;
    for (auto c_it = disjuncts.begin(); c_it != disjuncts.end(); ++c_it)
    {
      state_formulas::state_formula disjunct = (*this).apply(*c_it);
      //mCRL2log(log::verbose) << " * disjunct: " << pp(disjunct) << std::endl;
      std::vector<state_formulas::state_formula> parts = accessors::split_conjuncts(disjunct);
      std::set<data::variable> substitution_vars;
      bool apply = false;
      data::mutable_map_substitution< std::map< data::variable, data::data_expression > > substitutions;
      for(auto v_it = x_variables.begin(); v_it != x_variables.end(); ++v_it)
      {
        bool is_point = false;
        data::data_expression expr = data::sort_bool::false_();
        for(auto p_it = parts.begin(); p_it != parts.end(); ++p_it)
        {
          auto p = accessors::find_point((*p_it), (*v_it), /* negative = */ false);
          is_point = p.first;
          if (is_point)
          {
            expr = p.second;
            //mCRL2log(log::debug) << "Substitute " << pp(*v_it) << " with expression: " << pp(expr) << std::endl;
            //FIXME: check for inconsistent points (e.g., forall x. x=3 && x=5)
            substitutions[(*v_it)] = expr;
            substitution_vars.insert(*v_it);
            apply = true;
          }
        }
      }
      // substitution [e/v]
      //mCRL2log(log::debug) << "disjunct: " << pp(disjunct) << std::endl;
      if (apply)
      {
        disjunct = state_formulas::replace_variables(disjunct, substitutions);
      }
      //mCRL2log(log::debug) << "rewriting..." << std::endl;
      disjunct = rewrite_formula(disjunct);
      //mCRL2log(log::debug) << "rewritten disjunct: " << pp(disjunct) << std::endl;
      // remove substituted variables from the set of declared variables
      std::set<data::variable> decl_vars;
      std::set_difference(x_variables.begin(), x_variables.end(),
          substitution_vars.begin(), substitution_vars.end(),
          std::inserter(decl_vars, decl_vars.end()));

      // remove unused variables (e.g., exists x. d => 3 ---> d => 3)
      std::set<data::variable> used_vars = state_formulas::find_free_variables(disjunct);
      std::set<data::variable> vars;
      std::set_intersection(decl_vars.begin(), decl_vars.end(),
          used_vars.begin(), used_vars.end(),
          std::inserter(vars, vars.end()));
      data::variable_list var_list(vars.begin(), vars.end());

      disjunct = optimized::exists(var_list, disjunct);
      if (tr::is_true(disjunct))
      {
        //if (mCRL2logEnabled(log::verbose) && x != free_body)
        //{
        //  mCRL2log(log::verbose) << "exists: " << pp(x) << std::endl;
        //  mCRL2log(log::verbose) << "  ----> " << pp(free_body) << std::endl;
        //}
        return free_body; // && tr::true_();
      }
      result_disjuncts.push_back(disjunct);
      //if (mCRL2logEnabled(log::debug) && (*c_it) != disjunct)
      //{
      //  mCRL2log(log::debug) << "  * disjunct: " << pp(*c_it) << std::endl;
      //  mCRL2log(log::debug) << "  *     ----> " << pp(disjunct) << std::endl;
      //}
    }
    (*this).leave(x);
    state_formulas::state_formula result =
        optimized::and_(
            free_body,
            optimized::join_or(result_disjuncts.begin(), result_disjuncts.end()));

    //if (mCRL2logEnabled(log::verbose) && x != result)
    //{
    //  mCRL2log(log::verbose) << "exists: " << pp(x) << std::endl;
    //  mCRL2log(log::verbose) << "  ----> " << pp(result) << std::endl;
    //}
    return result;
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
void rewrite_quantifiers(T& x, data::data_specification& data, typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  rewrite_quantifiers_builder f(data);
  f.apply(x);
}

/// \brief
/// \param x an object containing state formulas
template <typename T>
T rewrite_quantifiers(const T& x, data::data_specification& data, typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  rewrite_quantifiers_builder f(data);
  T result = f.apply(x);
  //mCRL2log(log::verbose) << std::endl << "rewrite_quantifiers: " << std::endl <<
  //                          "formula: "<< pp(x) << std::endl;
  //mCRL2log(log::verbose) << "  ----> " << pp(result) << std::endl;
  mCRL2log(log::verbose) << "(" << f.report_cache_hits() << ")" << std::endl;
  return result;
}


} // namespace state_formulas

} // namespace mcrl2


#endif /* MCRL2_MODAL_FORMULAS_REWRITE_QUANTIFIERS_H_ */
