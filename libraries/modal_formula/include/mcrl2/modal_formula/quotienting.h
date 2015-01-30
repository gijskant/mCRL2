// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/quotienting.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULAS_QUOTIENTING_H_
#define MCRL2_MODAL_FORMULAS_QUOTIENTING_H_

#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include "boost/static_assert.hpp"
#include "mcrl2/data/bool.h"
#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/synchronization_vector.h"
#include "mcrl2/modal_formula/builder.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/quotienting_utility.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/modal_formula/state_formula_traits.h"
#include "mcrl2/modal_formula/state_formulas_optimized.h"
#include "mcrl2/modal_formula/detail/action_formula_in.h"

#include "mcrl2/process/label_generator.h"
#include "mcrl2/utilities/detail/join.h"


namespace mcrl2 {

namespace state_formulas {

/// \cond INTERNAL_DOCS

template <typename T>
void quotient(T& x, const lps::specification& spec, const lps::synchronization_vector& v, size_t i, quotienting_options options = quotienting_options(),
    typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

template <typename T>
T quotient(const T& x, const lps::specification& spec, const lps::synchronization_vector& v, size_t i, quotienting_options options = quotienting_options(),
    typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type* = 0);

// \brief Visitor for quotienting a state formula.
struct quotient_builder: public state_formula_builder<quotient_builder>
{
  typedef state_formula_builder<quotient_builder> super;
  typedef core::term_traits_optimized<state_formulas::state_formula> optimized;
  typedef core::term_traits<state_formulas::state_formula> tr;
  using super::enter;
  using super::leave;
  using super::update;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/builder_msvc.inc.h"
#endif

  lps::specification spec;
  lps::synchronization_vector v;
  size_t i; // index of component to be quotiented out
  data::set_identifier_generator id_generator;
  process::label_generator m_label_generator;
  quotienting_options m_options;
  std::stack<core::identifier_string> fixpoint_variables;
  std::map<core::identifier_string, std::map<data::variable, data::variable> > generated_variables;
  std::map<state_formulas::state_formula,
    std::map< std::map<data::variable, data::data_expression>,
      state_formulas::state_formula > > rewrite_cache;
  std::map<core::identifier_string, lps::synchronization_vector_type> vector_map;

  bool equation_system_mode;

  size_t n_cache_calls;
  size_t n_cache_hits;

  state_formulas::state_formula replace_variables(const state_formulas::state_formula& x, const std::map<data::variable, data::data_expression>& s)
  {
    data::mutable_map_substitution< std::map< data::variable, data::data_expression > > sigma(s);
    if (tr::is_true(x) || tr::is_false(x))
    {
      return x;
    }
    else
    {
      if (m_options.use_rewrite_cache)
      {
        n_cache_calls++;
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
            n_cache_hits++;
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

  /// \brief Generates a substitution function for the rewriter.
  /// \param v A vector of data variables
  /// \param e A vector of data expressions
  /// \return A substitution function.
  std::map< data::variable, data::data_expression >  make_substitution(data::variable_vector v, data::data_expression_vector e)
  {
    assert(v.size() == e.size());
    std::map< data::variable, data::data_expression > sigma;
    data::variable_vector::iterator i = v.begin();
    data::data_expression_vector::iterator j = e.begin();
    for (; i != v.end(); ++i, ++j)
    {
      sigma[*i] = *j;
    }
    return sigma;
  }

  /// \brief Generates a substitution function for the rewriter.
  /// \param a A list of assignments
  /// \return A substitution function.
  std::map< data::variable, data::data_expression >  make_substitution(data::assignment_list a)
  {
    mCRL2log(log::debug) << "    make_substitution:" << std::endl;
    std::map< data::variable, data::data_expression > sigma;
    for (auto a_it = a.begin(); a_it != a.end(); ++a_it)
    {
      sigma[(*a_it).lhs()] = (*a_it).rhs();
      mCRL2log(log::debug) << "      " << pp((*a_it).lhs()) << " --> " << pp((*a_it).rhs()) << std::endl;
    }
    return sigma;
  }

private:
  void init()
  {
    mCRL2log(log::verbose) << "use rewrite cache: " << (m_options.use_rewrite_cache ? "true":"false") << std::endl;
    mCRL2log(log::verbose) << "use vector map: " << (m_options.use_vector_map ? "true":"false") << std::endl;
    mCRL2log(log::debug) << "quotient_builder: i = " << i << std::endl;
    std::set<core::identifier_string> ids;
    ids = lps::find_identifiers(spec);
    id_generator.add_identifiers(ids);
    ids = data::find_identifiers(spec.data().constructors());
    id_generator.add_identifiers(ids);
    ids = data::find_identifiers(spec.data().mappings());
    id_generator.add_identifiers(ids);
    std::set<process::action_label> labels = lps::find_action_labels(spec);
    for (auto l_it = labels.begin(); l_it != labels.end(); ++l_it)
    {
      m_label_generator.add_identifier((*l_it).name());
    }
    if (m_options.use_vector_map)
    {
      for (auto e: v.vector())
      {
        core::identifier_string a = e.second.name();
        auto it = vector_map.find(a);
        if (it == vector_map.end())
        {
          lps::synchronization_vector_type vector;
          vector.push_back(e);
          vector_map[a] = vector;
        }
        else
        {
          it->second.push_back(e);
        }
      }
      /*mCRL2log(log::verbose) << "vector_map:" << std::endl;
      for(auto it1 = vector_map.begin(); it1 != vector_map.end(); ++it1)
      {
        mCRL2log(log::verbose) << (it1->first) << ": " << std::endl;
        for(auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
          mCRL2log(log::verbose) << " - " << pp(*it2) << std::endl;
        }
      }*/
    }
  }
public:

  quotient_builder(const lps::specification& spec_, const lps::synchronization_vector& v_, int i_,
      process::label_generator label_generator,
      quotienting_options options = quotienting_options())
    : spec(spec_),
      v(v_),
      i(i_),
      m_label_generator(label_generator),
      m_options(options),
      equation_system_mode(false),
      n_cache_calls(0),
      n_cache_hits(0)
  {
    init();
  }

  quotient_builder(const lps::specification& spec_, const lps::synchronization_vector& v_, int i_,
      quotienting_options options = quotienting_options())
    : spec(spec_),
      v(v_),
      i(i_),
      m_options(options),
      equation_system_mode(false),
      n_cache_calls(0),
      n_cache_hits(0)
  {
    init();
  }

  void set_equation_system_mode(bool mode)
  {
    equation_system_mode = mode;
  }

  std::string report_cache_hits()
  {
    std::stringstream s;
    s <<  "quotienting cache: " << n_cache_hits << " / " << n_cache_calls << " hits (" <<
         ((n_cache_calls==0) ? 100 : (((float)n_cache_hits/(float)n_cache_calls)*100)) << " %)";
    return s.str();
  }

  process::label_generator& label_generator()
  {
    return m_label_generator;
  }

  void add_identifiers(const modal_equation_system& x)
  {
    for (auto e : x.equations())
    {
      add_identifiers(e);
    }
    add_identifiers(x.initial_state());
  }

  void add_identifiers(const modal_equation& x)
  {
    id_generator.add_identifier(x.name());
    m_label_generator.add_identifier(x.name());
    for (auto p : x.parameters())
    {
      id_generator.add_identifier(p.name());
      m_label_generator.add_identifier(p.name());
    }
    add_identifiers(x.formula());
  }

  void add_identifiers(const state_formulas::state_formula& x)
  {
    std::set<core::identifier_string> ids = state_formulas::find_identifiers(x);
    id_generator.add_identifiers(ids);
    for (auto id : ids)
    {
      m_label_generator.add_identifier(id);
    }
  }

  state_formulas::state_formula data_to_state(const data::data_expression& x)
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

  template <typename QuantifyOperator, typename CombinatorOperator>
  inline state_formulas::state_formula quantify_combine(const data::variable_list& e, const state_formulas::state_formula& left, const state_formulas::state_formula& right,
      QuantifyOperator quantify, CombinatorOperator combine)
  {
    state_formulas::state_formula result =
        quantify(e,
              combine(
                  left,
                  right));
    return result;
  }

  template <typename QuantifyOperator, typename CombinatorOperator, typename ModalOperator, typename JoinOperator>
  state_formulas::state_formula modality(const regular_formulas::regular_formula& x_formula, const state_formulas::state_formula& x_operand,
      QuantifyOperator quantify, CombinatorOperator combine, ModalOperator modal, JoinOperator join, state_formulas::state_formula zero)
  {
    std::set<state_formulas::state_formula> parts;
    std::set<std::pair<state_formulas::state_formula, data::variable_list>> onlyactive_formula_signatures;
    std::set<std::pair<state_formulas::state_formula, std::pair<action_formulas::action_formula, process::action_label>>> inactive_formula_signatures;
    assert(action_formulas::is_action_formula(x_formula));
    const action_formulas::action_formula& alpha = atermpp::down_cast<const action_formulas::action_formula>(x_formula);
    mCRL2log(log::debug) << "  action formula: " << pp(alpha) << std::endl;
    mCRL2log(log::debug) << "  operand: " << pp(x_operand) << std::endl;
    mCRL2log(log::debug) << "  computing quotient formula for operand..." << std::endl;
    state_formulas::state_formula phi_quotient = (*this).apply(x_operand);
    mCRL2log(log::debug) << "  phi_quotient: " << pp(phi_quotient) << std::endl;

    //mCRL2log(log::debug) << "  computing subformulae for action formula..." << std::endl;
    lps::synchronization_vector_type vector;
    bool use_vector_map = false;
    core::identifier_string formula_action_name;
    if (m_options.use_vector_map && action_formulas::is_multi_action(alpha))
    {
      const action_formulas::multi_action& multi = atermpp::down_cast<const action_formulas::multi_action>(alpha);
      if (multi.actions().size()==1)
      {
        const process::action& act = multi.actions().front();
        formula_action_name = act.label().name();
        use_vector_map = true;
      }
    }
    if (use_vector_map)
    {
      auto it = vector_map.find(formula_action_name);
      if (it != vector_map.end())
      {
        vector = it->second;
      }
    }
    else
    {
      vector= v.vector();
    }
    for (auto v_it = vector.begin(); v_it != vector.end(); ++v_it) {
      std::string a_i = (*v_it).first[i];
      process::action_label a((*v_it).second);
      //mCRL2log(log::debug) << "  [ vector element: " << pp(a) << " ";
      //if (mCRL2logEnabled(log::debug))
      //{
      //  std::copy((*v_it).first.begin(), (*v_it).first.end(), std::ostream_iterator<std::string>(mcrl2::log::mcrl2_logger().get(log::debug), " "));
      //}
      //mCRL2log(log::debug) << " (" << a_i << ") ]" << std::endl;
      if (a_i == lps::inactive_label) {
        // case P_i is inactive
        //mCRL2log(log::verbose) << "  * P_i is inactive" << std::endl;
        auto formula_signature = std::make_pair(phi_quotient, std::make_pair(alpha, a));
        if (inactive_formula_signatures.find(formula_signature) != inactive_formula_signatures.end())
        {
          // skip
          //mCRL2log(log::verbose) << "  - skip formula for: phi quotient = " << pp(phi_quotient) << ", alpha = " << pp(alpha) << ", action label = " << pp(a) << " (" << pp(a.sorts()) << ")" << std::endl;
        }
        else
        {
          inactive_formula_signatures.insert(formula_signature);
          //mCRL2log(log::verbose) << "  - add formula for: phi quotient = " << pp(phi_quotient) << ", alpha = " << pp(alpha) << ", action label = " << pp(a) << " (" << pp(a.sorts()) << ")" << std::endl;
          std::vector<data::variable> l;
          for (auto s = a.sorts().begin(); s != a.sorts().end(); ++s)
          {
            l.push_back(data::variable(id_generator(std::string("d")), (*s)));
          }
          process::action a_action(a, data::data_expression_list(l.begin(), l.end()));
          process::action_list a_action_list;
          a_action_list.push_front(a_action);
          action_formulas::multi_action a_multi_action(a_action_list);

          state_formulas::state_formula subformula =
              quantify_combine(
                  data::variable_list(l.begin(), l.end()),
                  data_to_state(action_formulas::detail::in(a_action, alpha, id_generator)),
                  modal(a_multi_action, phi_quotient),
                  quantify,
                  combine);
          //if (mCRL2logEnabled(log::debug))
          //{
          //  mCRL2log(log::debug) << "  - subformula: " << pp(subformula) << std::endl;
          //}
          if (subformula==zero)
          {
            return zero;
          }
          parts.insert(subformula);
        }
      }
      else
      {
        bool other_process_active = false;
        //mCRL2log(log::debug) << "  * P_i not inactive. Determining if other processes are active." << std::endl;
        for (size_t j = 0; j < v.length(); j++)
        {
          if (i != j && (*v_it).first[j] != lps::inactive_label)
          {
            //mCRL2log(log::debug) << "  (Process " << j << " is active with label " << (*v_it).first[j] << ")" << std::endl;
            other_process_active = true; break;
          }
        }
        if (!other_process_active)
        {
          // case P_i is the only active process
          //mCRL2log(log::verbose) << "  * P_i is only active" << std::endl;
          for(auto summand_it = spec.process().action_summands().begin(); summand_it != spec.process().action_summands().end(); ++ summand_it)
          {
            lps::action_summand summand = (*summand_it);
            lps::multi_action multi = summand.multi_action();
            process::action action;
            if (multi.actions().empty()) {
              action = process::action(process::action_label("tau", data::sort_expression_list()), data::data_expression_list());
            } else if (multi.actions().size() == 1) {
              action = multi.actions().front();
            } else {
              throw std::runtime_error("Quotienting is only defined for singleton actions: " + pp(multi));
            }
            //mCRL2log(log::verbose) << "  # summand action = " << pp(action) << std::endl;

            std::string label_name = action.label().name();
            if (a_i == label_name && a.sorts() == action.label().sorts())
            {
              // substitution [g_{1,i}/d_1]
              mCRL2log(log::debug) << "    make_substitution:" << std::endl;
              std::map< data::variable, data::data_expression > substitution;
              for (auto a: summand.assignments())
              {
                if (!equation_system_mode)
                {
                  substitution[generated_variables[fixpoint_variables.top()][a.lhs()]] = a.rhs();
                  mCRL2log(log::debug) << "      " << pp(generated_variables[fixpoint_variables.top()][a.lhs()]) << " --> " << pp(a.rhs()) << std::endl;
                }
                substitution[a.lhs()] = a.rhs();
                mCRL2log(log::debug) << "      " << pp(a.lhs()) << " --> " << pp(a.rhs()) << std::endl;
              }
              state_formulas::state_formula phi_quotient_subst = replace_variables(phi_quotient, substitution);

              mCRL2log(log::debug) << "modality: " << std::endl <<
                  "before substitution: " << phi_quotient << std::endl <<
                  "            -------> " << phi_quotient_subst << std::endl << std::endl;

              data::data_expression c = summand.condition();
              process::action a_action(a, action.arguments());
              state_formulas::state_formula body =
                  combine(
                      optimized::and_(
                          data_to_state(action_formulas::detail::in(a_action, alpha, id_generator)),
                          data_to_state(c)),
                      phi_quotient_subst);

              // FIXME: check if no earlier subformula with the same body has been added.
              // To prevent generation of formulas like:
              // forall a1:A, b1:B . (c -> X(d, e))[d:=a1, e:=b1]
              // && forall a2:A, b2:B . (c -> X(d, e))[d:=a2, e:=b2]
              // && ...
              auto formula_signature = std::make_pair(body, summand.summation_variables());
              if (onlyactive_formula_signatures.find(formula_signature) != onlyactive_formula_signatures.end())
              {
                // skip
                //if (mCRL2logEnabled(log::verbose))
                //{
                //  mCRL2log(log::verbose) << "  - skip body: " << pp(body) << " for summand (with variables: " << pp(summand.summation_variables()) << ")"<< std::endl;
                //}
              }
              else
              {
                onlyactive_formula_signatures.insert(formula_signature);
                //mCRL2log(log::verbose) << "  - add body: " << pp(body) << " for summand (with variables: " << pp(summand.summation_variables()) << ")"<< std::endl;
                // generate quantifier variables and build substitution
                std::vector<data::variable> new_summation_variables;
                std::vector<data::variable> match_v;
                std::vector<data::data_expression> subst_v;
                for (auto s_v_it = summand.summation_variables().begin(); s_v_it != summand.summation_variables().end(); ++s_v_it)
                {
                  data::variable summation_variable = *s_v_it;
                  data::variable new_variable = data::variable(id_generator(summation_variable.name()), summation_variable.sort());
                  new_summation_variables.push_back(new_variable);
                  match_v.push_back(summation_variable);
                  subst_v.push_back(new_variable);
                }
                std::map< data::variable, data::data_expression > summation_variable_subst =
                                  make_substitution(match_v, subst_v);
                data::variable_list e(new_summation_variables.begin(), new_summation_variables.end());

                state_formulas::state_formula body_subst = replace_variables(body, summation_variable_subst);
                //if (mCRL2logEnabled(log::debug))
                //{
                //  mCRL2log(log::debug) << "body: " << pp(body) << std::endl << "body_subst: " << pp(body_subst) << std::endl;
                //}
                state_formulas::state_formula subformula =
                    quantify(e, body_subst);

                if (subformula==zero)
                {
                  return zero;
                }
                parts.insert(subformula);
                // FIXME: mark subformula for unfolding
                //if (mCRL2logEnabled(log::debug))
                //{
                //  mCRL2log(log::debug) << "  - subformula: " << pp(subformula) << std::endl;
                //}
              }
            }
          }
        }
        else
        {
          // case P_i synchronises with other processes
          //mCRL2log(log::debug) << "  * P_i synchronises" << std::endl;
          for(auto summand_it = spec.process().action_summands().begin(); summand_it != spec.process().action_summands().end(); ++ summand_it)
          {
            lps::action_summand summand = (*summand_it);
            lps::multi_action multi = summand.multi_action();
            process::action action;
            if (multi.actions().empty()) {
              action = process::action(process::action_label("tau", data::sort_expression_list()), data::data_expression_list());
            } else if (multi.actions().size() == 1) {
              action = multi.actions().front();
            } else {
              throw std::runtime_error("Quotienting is only defined for singleton actions: " + pp(multi));
            }
            //mCRL2log(log::debug) << "  # summand action = " << pp(action) << std::endl;

            std::string label_name = action.label().name();
            if (a_i == label_name && a.sorts() == action.label().sorts())
            {
              // substitution [g_{1,i}/d_1]
              std::map< data::variable, data::data_expression > substitution =
                  make_substitution(summand.assignments());
              state_formula phi_quotient_subst = replace_variables(phi_quotient, substitution);

              data::variable_list e = summand.summation_variables();
              data::data_expression c = summand.condition();
              process::action_label rho_label = m_label_generator.fresh_action_label((*v_it).first, a);
              process::action_list rho_actions;
              rho_actions.push_front(process::action(rho_label, action.arguments()));
              action_formulas::multi_action rho_action(rho_actions);
              action_formulas::action_formula rho_formula(rho_action);
              regular_formulas::regular_formula rho(rho_formula);
              process::action a_action(a, action.arguments());
              state_formulas::state_formula subformula =
                  quantify_combine(
                      e,
                      optimized::and_(
                          data_to_state(action_formulas::detail::in(a_action, alpha, id_generator)),
                          data_to_state(c)),
                      modal(rho, phi_quotient_subst),
                      quantify,
                      combine);
              if (subformula==zero)
              {
                return zero;
              }
              parts.insert(subformula);
              //if (mCRL2logEnabled(log::debug))
              //{
              //  mCRL2log(log::debug) << "  - subformula: " << pp(subformula) << std::endl;
              //}
            }
          }
        }
      }
    }
    //if (skipped_subformulas > 0)
    //{
    //  mCRL2log(log::debug) << "Skipped " << skipped_subformulas << " redundant subformulas" << std::endl;
    //}
    state_formulas::state_formula result = join(parts.begin(), parts.end());
    if (is_true(result) || is_false(result))
    {
      // skip unfolding // FIXME
    }
    return result;
  }

  static inline state_formulas::state_formula must(const regular_formulas::regular_formula& formula, const state_formula& operand)
  {
    return state_formulas::must(formula, operand);
  }

  static inline state_formulas::state_formula may(const regular_formulas::regular_formula& formula, const state_formula& operand)
  {
    return state_formulas::may(formula, operand);
  }

  static inline
  state_formulas::state_formula imp(const state_formulas::state_formula& p, const state_formulas::state_formula& q)
  {
    if (tr::is_data(p)) {
      data::data_expression l(p);
      return optimized::or_(state_formulas::state_formula(data::sort_bool::not_(l)), q);
    } else {
      return optimized::or_(tr::not_(p), q);
    }
  }

  state_formulas::state_formula apply(const state_formulas::must& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "must [" << pp(x.formula()) << "] " << pp(x.operand()) << ":" << std::endl;
    state_formulas::state_formula result = modality(x.formula(), x.operand(),
        optimized::forall, imp, &must, optimized::join_and<std::set<state_formulas::state_formula>::iterator >, tr::false_());
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::may& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "may <" << pp(x.formula()) << "> " << pp(x.operand()) << ":" << std::endl;
    state_formulas::state_formula result = modality(x.formula(), x.operand(),
        optimized::exists, optimized::and_, &may, optimized::join_or<std::set<state_formulas::state_formula>::iterator >, tr::true_());
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::variable& x)
  {
    (*this).enter(x);
    std::vector<data::data_expression> arguments(x.arguments().begin(), x.arguments().end());
    arguments.insert(arguments.end(), spec.process().process_parameters().begin(), spec.process().process_parameters().end());
    (*this).leave(x);
    state_formulas::state_formula result = state_formulas::variable(x.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  template <typename FixpointOperator>
  inline state_formulas::state_formula fixpoint(const core::identifier_string& x_name,
      const data::assignment_list& x_assignments,
      const state_formulas::state_formula x_operand,
      FixpointOperator sigma
      )
  {
    if (equation_system_mode)
    {
      throw std::runtime_error("fixpoint formula in equation system mode.");
    }
    //mCRL2log(log::debug) << "fixpoint formula: name = " << x_name << std::endl;
    //mCRL2log(log::debug) << "  assignments = " << pp(x_assignments) << std::endl;
    data::assignment_vector assignments(x_assignments.begin(), x_assignments.end());
    // create list of substitutions
    data::variable_vector substitution_variables;
    data::data_expression_vector substitutions;
    // add assignments for process parameters.
    for(auto it = spec.initial_process().assignments().begin(); it != spec.initial_process().assignments().end(); ++it)
    {
      data::assignment a = *it;
      data::variable d = a.lhs();
      // keep a stack of fixpoint variables, use it in generating the name
      core::identifier_string name = d.name();
      if (!fixpoint_variables.empty()) {
        std::ostringstream s;
        s << name << "_" << fixpoint_variables.top();
        name = core::identifier_string(s.str());
      }
      name = id_generator(name);
      data::variable d_x(name, d.sort());
      generated_variables[x_name][d] = d_x;
      data::data_expression e_x;
      if (fixpoint_variables.empty()) {
        // if Y is the outermost fixpoint variable, use the assigment from the LPS;
        e_x = a.rhs();
      } else {
        // otherwise, instead of rhs, use d_1^X (see paper).
        e_x = generated_variables[fixpoint_variables.top()][d];
      }
      data::assignment a_x(d_x, e_x);
      assignments.push_back(a_x);
      substitution_variables.push_back(d);
      substitutions.push_back(d_x);
      //mCRL2log(log::debug) << "  - assignment: " << pp(a_x) << ", variable: " << pp(d) << ", substitution: " << pp(d_x) << std::endl;
    }
    fixpoint_variables.push(x_name); // push variable Y to stack
    std::map< data::variable, data::data_expression > substitution =
        make_substitution(substitution_variables, substitutions);
    //mCRL2log(log::debug) << "  new assignments = " << pp(data::assignment_list(assignments.begin(), assignments.end())) << std::endl;
    //mCRL2log(log::debug) << "  operand = " << pp(x_operand) << std::endl;
    //mCRL2log(log::debug) << "  computing quotient for operand." << std::endl;
    state_formulas::state_formula phi_quotient =  (*this).apply(x_operand); // phi //_Y P_1
    //mCRL2log(log::debug) << "  phi_quotient = " << pp(phi_quotient) << std::endl;
    phi_quotient = replace_variables(phi_quotient, substitution);  // [d_1^Y / d_1]
    // todo: substitute data variables
    state_formulas::state_formula result = sigma(x_name, data::assignment_list(assignments.begin(), assignments.end()), phi_quotient);
    fixpoint_variables.pop(); // pop variable Y
    return result;
  }

  static inline state_formulas::state_formula nu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand)
  {
    return state_formulas::nu(name, assignments, operand);
  }

  static inline state_formulas::state_formula mu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand)
  {
    return state_formulas::mu(name, assignments, operand);
  }

  state_formulas::state_formula apply(const state_formulas::nu& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "nu " << x.name() << ": " << pp(x) << std::endl;
    state_formulas::state_formula result = fixpoint(x.name(), x.assignments(), x.operand(), &nu);
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::mu& x)
  {
    (*this).enter(x);
    //mCRL2log(log::debug) << "mu " << x.name() << ": " << pp(x) << std::endl;
    state_formulas::state_formula result = fixpoint(x.name(), x.assignments(), x.operand(), &mu);
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
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
    //mCRL2log(log::debug) << "and: " << x << std::endl;
    state_formulas::state_formula result = optimized::and_((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::or_& x)
  {
    (*this).enter(x);
    state_formulas::state_formula result = optimized::or_((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::imp& x)
  {
    (*this).enter(x);
    state_formulas::state_formula result = optimized::imp((*this).apply(x.left()), (*this).apply(x.right()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::forall& x)
  {
    (*this).enter(x);
    state_formulas::state_formula result = optimized::forall(x.variables(), (*this).apply(x.body()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  state_formulas::state_formula apply(const state_formulas::exists& x)
  {
    (*this).enter(x);
    state_formulas::state_formula result = optimized::exists(x.variables(), (*this).apply(x.body()));
    (*this).leave(x);
    //mCRL2log(log::debug) << "result: " << pp(result) << std::endl;
    return result;
  }

  void update(state_formulas::modal_equation& x)
  {
    mCRL2log(log::verbose) << "Quotienting equation " << x.name() << std::endl;
    std::vector<data::variable> parameters(x.parameters().begin(), x.parameters().end());
    // add process parameters.
    std::map< data::variable, data::data_expression > substitution;
    for(auto a : spec.initial_process().assignments())
    {
      core::identifier_string name = a.lhs().name();
      {
        std::ostringstream s;
        s << name << "_" << x.name();
        name = core::identifier_string(s.str());
      }
      name = id_generator(name);
      data::variable d_x(name, a.lhs().sort());
      parameters.push_back(d_x);
      substitution[a.lhs()] = d_x;
    }
    x.set_parameters(data::variable_list(parameters.begin(), parameters.end()));
    state_formulas::state_formula phi_quotient =  (*this).apply(x.formula()); // phi //_Y P_1
    phi_quotient = replace_variables(phi_quotient, substitution);  // [d_1^Y / d_1]
    x.set_formula(phi_quotient);
    //FIXME: apply unfold_recursion on marked parts
  }

  void update(state_formulas::modal_equation_system& x)
  {
    set_equation_system_mode(true);
    std::vector<modal_equation> result_equations;
    for(auto e : x.equations())
    {
      (*this).update(e);
      result_equations.push_back(e);
    }
    variable initial_variable = x.initial_state();
    std::vector<data::data_expression> arguments(initial_variable.arguments().begin(), initial_variable.arguments().end());
    for(auto a : spec.initial_process().assignments())
    {
      arguments.push_back(a.rhs());
    }
    initial_variable = variable(initial_variable.name(), data::data_expression_list(arguments.begin(), arguments.end()));
    x.set_initial_state(initial_variable);
    x.set_equations(result_equations);
  }

};
/// \endcond

/// \brief
/// \param x an object containing state formulas
template <typename T>
void quotient(T& x, const lps::specification& spec, const lps::synchronization_vector& v, size_t i, quotienting_options options,
    typename std::enable_if<!std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  quotient_builder f(spec, v, i, options);
  f.add_identifiers(x);
  f.update(x);
}

/// \brief
/// \param x an object containing state formulas
template <typename T>
T quotient(const T& x, const lps::specification& spec, const lps::synchronization_vector& v, size_t i, quotienting_options options,
    typename std::enable_if<std::is_base_of<atermpp::aterm, T>::value>::type*)
{
  quotient_builder f(spec, v, i, options);
  f.add_identifiers(x);
  T result = f.apply(x);
  mCRL2log(log::verbose) << "(" << f.report_cache_hits() << ")" << std::endl;
  return result;
}

} // namespace state_formulas

} // namespace mcrl2

#endif /* MCRL2_MODAL_FORMULAS_QUOTIENTING_H_ */
