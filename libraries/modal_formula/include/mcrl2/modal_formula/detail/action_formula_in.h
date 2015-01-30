// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/detail/action_formula_in.h
/// \brief add your file description here.
/// Based on the code of pbes/detail/lps2pbes_sat.h

#ifndef MCRL2_MODAL_FORMULA_DETAIL_ACTION_FORMULA_IN_H
#define MCRL2_MODAL_FORMULA_DETAIL_ACTION_FORMULA_IN_H

#include <string>

#include "mcrl2/data/standard_utility.h"
#include "mcrl2/data/detail/data_utility.h"
#include "mcrl2/data/detail/find.h"
#include "mcrl2/modal_formula/replace.h"
#include "mcrl2/modal_formula/traverser.h"
//#include "mcrl2/pbes/detail/lps2pbes_utility.h"
#include "mcrl2/data/optimized_boolean_operators.h"


namespace mcrl2 {

namespace action_formulas {

namespace detail {

/// \brief Generates a substitution that assigns fresh variables to the given sequence of variables.
/// The identifier generator is used to assign names to the fresh variables.
/// Caveat: the implementation is very inefficient.
/// \param update_context If true, then generated names are added to the context
inline
data::mutable_map_substitution<> make_fresh_variables(const data::variable_list& variables, data::set_identifier_generator& id_generator, bool add_to_context = true)
{
  data::mutable_map_substitution<> result;
  for (auto i = variables.begin(); i != variables.end(); ++i)
  {
    core::identifier_string name =  id_generator(std::string(i->name()));
    result[*i] = data::variable(name, i->sort());
    if (!add_to_context)
    {
      id_generator.remove_identifier(name);
    }
  }
  return result;
}

/// \brief
data::data_expression in(const process::action& a,
                    const action_formulas::action_formula& x,
                    data::set_identifier_generator& id_generator
                   );

template <typename Derived>
struct in_traverser: public action_formulas::action_formula_traverser<Derived>
{
  typedef action_formulas::action_formula_traverser<Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif

  const process::action& a;
  data::set_identifier_generator& id_generator;
  std::vector<data::data_expression> result_stack;

  in_traverser(const process::action& a_, data::set_identifier_generator& id_generator_)
    : a(a_), id_generator(id_generator_)
  {}

  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  void push(const data::data_expression& x)
  {
    result_stack.push_back(x);
  }

  const data::data_expression& top() const
  {
    return result_stack.back();
  }

  data::data_expression pop()
  {
    data::data_expression result = top();
    result_stack.pop_back();
    return result;
  }

  void leave(const data::data_expression& x)
  {
    push(x);
  }

  bool is_tau(process::action_label lps_action)
  {
    std::string lps_action_name(lps_action.name());
    if (lps_action_name.size() >= 3 && lps_action_name.substr(0,3) == "tau")
    {
      std::string nr = lps_action_name.substr(3, lps_action_name.size()-1);
      if (std::atoi(nr.c_str())==0)
      {
        return false;
      }
      return true;
    }
    return false;
  }

  void leave(const action_formulas::multi_action& x)
  {
    process::action action;
    bool match_tau = false;
    if (x.actions().size() == 0)
    {
      action = process::action(process::action_label("tau", data::sort_expression_list()), data::data_expression_list());
      if (is_tau(a.label()))
      {
        match_tau = true;
      }
    }
    else if (x.actions().size() > 1)
    {
      throw std::runtime_error("Multi-actions not supported.");
    }
    else
    {
      action = x.actions().front();
    }
    if (match_tau)
    {
      push(data::sort_bool::true_());
    }
    else if (a.label()==action.label())
    {
      std::vector<data::data_expression> subformulae;
      assert(a.arguments().size()==action.arguments().size());
      auto it1 = a.arguments().begin();
      auto it2 = action.arguments().begin();
      for(; it1 != a.arguments().end(); ++it1, ++it2)
      {
        subformulae.push_back(data::lazy::equal_to((*it1), (*it2)));
      }
      push(data::lazy::join_and(subformulae.begin(), subformulae.end()));
    }
    else
    {
      push(data::sort_bool::false_());
    }
  }

  void leave(const action_formulas::true_&)
  {
    push(data::sort_bool::true_());
  }

  void leave(const action_formulas::false_&)
  {
    push(data::sort_bool::false_());
  }

  void apply(const action_formulas::not_& x)
  {
    push(data::lazy::not_(in(a, x.operand(), id_generator)));
  }

  void leave(const action_formulas::and_&)
  {
    data::data_expression right = pop();
    data::data_expression left = pop();
    push(data::lazy::and_(left, right));
  }

  void leave(const action_formulas::or_&)
  {
    data::data_expression right = pop();
    data::data_expression left = pop();
    push(data::optimized_or(left, right));
  }

  void leave(const action_formulas::imp&)
  {
    data::data_expression right = pop();
    data::data_expression left = pop();
    push(data::lazy::implies(left, right));
  }

  void apply(const action_formulas::forall& x)
  {
    data::mutable_map_substitution<> sigma_x = make_fresh_variables(x.variables(), id_generator, false);
    std::set<data::variable> sigma_x_variables = data::substitution_variables(sigma_x);
    action_formulas::action_formula alpha = x.body();
    data::variable_list y = data::replace_variables(x.variables(), sigma_x);
    push(data::optimized_forall(y, in(a, action_formulas::replace_variables_capture_avoiding(alpha, sigma_x, sigma_x_variables), id_generator)));
  }

  void apply(const action_formulas::exists& x)
  {
    data::mutable_map_substitution<> sigma_x = make_fresh_variables(x.variables(), id_generator, false);
    std::set<data::variable> sigma_x_variables = data::substitution_variables(sigma_x);
    action_formulas::action_formula alpha = x.body();
    data::variable_list y = data::replace_variables(x.variables(), sigma_x);
    push(data::optimized_exists(y, in(a, action_formulas::replace_variables_capture_avoiding(alpha, sigma_x, sigma_x_variables), id_generator)));
  }

};

template <template <class> class Traverser>
struct apply_in_traverser: public Traverser<apply_in_traverser<Traverser> >
{
  typedef Traverser<apply_in_traverser<Traverser> > super;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::top;

  apply_in_traverser(const process::action& a, data::set_identifier_generator& id_generator)
    : super(a, id_generator)
  {}

#ifdef BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif
};

} // namespace detail

} // namespace action_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_DETAIL_ACTION_FORMULA_IN_H
