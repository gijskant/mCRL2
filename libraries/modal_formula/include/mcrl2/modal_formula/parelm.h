// Author(s): Gijs Kant
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/parelm.h
/// \brief The parelm algorithm.

#ifndef MCRL2_MODAL_FORMULA_PARELM_H
#define MCRL2_MODAL_FORMULA_PARELM_H

#include "mcrl2/utilities/reachable_nodes.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/utilities/detail/iota.h"
#include "mcrl2/modal_formula/algorithms.h"
#include "mcrl2/modal_formula/modal_equation_system.h"
#include "mcrl2/modal_formula/traverser.h"
#include <iterator>
#include <boost/graph/adjacency_list.hpp>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <map>
#include <set>
#include <vector>

namespace mcrl2
{

namespace state_formulas
{

template <typename Derived>
struct parelm_traverser: public state_formulas::state_formula_traverser<Derived>
{
  typedef state_formulas::state_formula_traverser<Derived> super;
  using super::enter;
  using super::leave;
  using super::apply;

  /// \brief The graph type of the dependency graph
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> graph;
  /// \brief The vertex type of the dependency graph
  typedef boost::graph_traits<graph>::vertex_descriptor vertex_descriptor;
  /// \brief The edge type of the dependency graph
  typedef boost::graph_traits<graph>::edge_descriptor edge_descriptor;


#if BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif

protected:

  size_t offset;
  std::vector<core::identifier_string> propvars;
  std::map<core::identifier_string, size_t> propvar_offset_start;
  std::map<core::identifier_string, size_t> propvar_offset_end;
  std::vector<data::variable> m_predicate_variables;
  std::stack<core::identifier_string> variable_stack;
  std::stack<std::set<data::variable> > bound_data_variable_stack;
  std::stack<std::map<data::variable, size_t> > variable_bindings_stack;
  const size_t N; // number of vertices in G
  graph G; // the dependency graph G
  std::set<size_t> m_markings;

public:

  parelm_traverser(size_t _N) :
    offset(0),
    N(_N),
    G(_N)
  {}

  const graph& dependency_graph()
  {
    return G;
  }

  const std::set<size_t>& markings()
  {
    return m_markings;
  }

  std::vector<data::variable>& predicate_variables()
  {
    return m_predicate_variables;
  }

  const std::vector<core::identifier_string>& variables()
  {
    return propvars;
  }

  std::map<core::identifier_string, size_t>& offset_start()
  {
    return propvar_offset_start;
  }

  std::map<core::identifier_string, size_t>& offset_end()
  {
    return propvar_offset_end;
  }


  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  static inline std::string print_variables(const std::set<data::variable>& variables)
  {
    std::stringstream str;
    str << "{";
    bool first = true;
    for(auto it = variables.begin(); it != variables.end(); ++it)
    {
      if (!first)
      {
        str << ", ";
      }
      str << data::pp(*it);
      first = false;
    }
    str << "}";
    return str.str();
  }

  inline std::string print_index(size_t i)
  {
    std::stringstream str;
    for(auto it = propvars.begin(); it != propvars.end(); ++it)
    {
      core::identifier_string name = (*it);
      size_t offset_end = propvar_offset_end[name];
      if (i < offset_end)
      {
        size_t offset_start = propvar_offset_start[name];
        str << "var: " << name << ", index: " << (i - offset_start) << ", param: " << pp(m_predicate_variables[i]);
        return str.str();
      }
    }
    throw runtime_error("print_index: could not find variable for index " + boost::lexical_cast<std::string>(i));
    return str.str();
  }

  void visit_free_variables(const std::set<data::variable>& x)
  {
    if (!variable_bindings_stack.empty())
    {
      std::set<data::variable> vars;
      if (bound_data_variable_stack.empty())
      {
        vars = x;
      }
      else
      {
        std::set_difference(x.begin(), x.end(),
            bound_data_variable_stack.top().begin(), bound_data_variable_stack.top().end(),
            std::inserter(vars, vars.end()));
      }
      for (auto v_it = vars.begin(); v_it != vars.end(); ++v_it)
      {
        data::variable var = (*v_it);
        //mCRL2log(log::debug) << "  - Looking up index for data variable " << pp(var) << std::endl;
        auto f_it = variable_bindings_stack.top().find(var);
        if (f_it == variable_bindings_stack.top().end())
        {
          //mCRL2log(log::debug) << "    not found." << std::endl;
        }
        else
        {
          size_t index = f_it->second;
          //mCRL2log(log::debug) << "    index = " << index << " (" << print_index(index) << ") [MARKED]" << std::endl;
          m_markings.insert(index);
        }
      }
    }
  }

  void apply(const data::data_expression& x)
  {
    // determine if state variables are 'used':
    // variables in x - bound variables.
    // then connect to variables of variable_stack.top()
    //mCRL2log(log::debug) << "Computing dependencies for data expression: " << pp(x) << std::endl;
    std::set<data::variable> free_variables = data::find_free_variables(x);
    //mCRL2log(log::debug) << "  free variables: ";
    //mCRL2log(log::debug) << print_variables(free_variables);
    visit_free_variables(free_variables);
  }

  void apply(const state_formulas::variable& x)
  {
    //mCRL2log(log::debug) << "Visit variable: " << pp(x) << std::endl;
    size_t local_offset = propvar_offset_start[x.name()];
    x.arguments();
    for (auto it = x.arguments().begin(); it != x.arguments().end(); ++it)
    {
      data::data_expression expr = (*it);
      if (!variable_bindings_stack.empty())
      {
        //mCRL2log(log::debug) << "  Computing dependencies for data expression: " << pp(expr)
        //    << " at offset " << local_offset << ", (variable " << core::pp(x.name()) << ")" << std::endl;
        std::set<data::variable> free_variables = data::find_free_variables(expr);
        std::set<data::variable> vars;
        if (bound_data_variable_stack.empty())
        {
          vars = free_variables;
        }
        else
        {
          std::set_difference(free_variables.begin(), free_variables.end(),
              bound_data_variable_stack.top().begin(), bound_data_variable_stack.top().end(),
              std::inserter(vars, vars.end()));
        }
        for (auto v_it = vars.begin(); v_it != vars.end(); ++v_it)
        {
          data::variable var = (*v_it);
          //mCRL2log(log::debug) << "  - Looking up index for data variable " << pp(var) << std::endl;
          auto f_it = variable_bindings_stack.top().find(var);
          if (f_it == variable_bindings_stack.top().end())
          {
            //mCRL2log(log::debug) << "    not found." << std::endl;
          }
          else
          {
            size_t index = f_it->second;
            //mCRL2log(log::debug) << "    index = " << index << " (" << print_index(index) << ")" << std::endl;
            //mCRL2log(log::debug) << "    => (" << print_index(index) << ") affects (" << print_index(local_offset) << ")" << std::endl;
            boost::add_edge(local_offset, index, G); // local_offset depends on index
          }
        }
      }
      local_offset++;
    }

  }

  void enter(const state_formulas::forall& x)
  {
    //mCRL2log(log::debug) << "Visit forall." << std::endl;
    std::set<data::variable> bound_variables;
    if (!bound_data_variable_stack.empty())
    {
      bound_variables = bound_data_variable_stack.top();
    }
    bound_variables.insert(x.variables().begin(), x.variables().end());
    bound_data_variable_stack.push(bound_variables);
  }

  void leave(const state_formulas::forall& x)
  {
    bound_data_variable_stack.pop();
  }

  void enter(const state_formulas::exists& x)
  {
    mCRL2log(log::debug) << "Visit exists." << std::endl;
    std::set<data::variable> bound_variables;
    if (!bound_data_variable_stack.empty())
    {
      bound_variables = bound_data_variable_stack.top();
    }
    bound_variables.insert(x.variables().begin(), x.variables().end());
    bound_data_variable_stack.push(bound_variables);
  }

  void leave(const state_formulas::exists& x)
  {
    bound_data_variable_stack.pop();
  }

  void push_variable(const core::identifier_string& name, const data::assignment_list& assignments)
  {
    //mCRL2log(log::debug) << "Visit fixpoint formula: " << core::pp(name) << std::endl;
    variable_stack.push(name);
    propvars.push_back(name);
    propvar_offset_start[name] = offset;
    propvar_offset_end[name] = offset + assignments.size();
    std::map<data::variable, size_t> variable_bindings;
    if (!variable_bindings_stack.empty())
    {
      variable_bindings = variable_bindings_stack.top();
    }
    //offset += assignments.size();
    for (auto it = assignments.begin(); it != assignments.end(); ++it)
    {
      data::assignment a = (*it);
      m_predicate_variables.push_back(a.lhs());
      variable_bindings[a.lhs()] = offset;
      if (!variable_bindings_stack.empty())
      {
        //mCRL2log(log::debug) << "Computing dependencies for data assignment: " << pp(a.lhs())
        //    << " at offset " << offset << ", (variable " << core::pp(name) << ")" << std::endl;
        std::set<data::variable> free_variables = data::find_free_variables(a.rhs());
        std::set<data::variable> vars;
        if (bound_data_variable_stack.empty())
        {
          vars = free_variables;
        }
        else
        {
          std::set_difference(free_variables.begin(), free_variables.end(),
              bound_data_variable_stack.top().begin(), bound_data_variable_stack.top().end(),
              std::inserter(vars, vars.end()));
        }
        for (auto v_it = vars.begin(); v_it != vars.end(); ++v_it)
        {
          data::variable var = (*v_it);
          //mCRL2log(log::debug) << "- Looking up index for data variable " << pp(var) << std::endl;
          auto f_it = variable_bindings_stack.top().find(var);
          if (f_it == variable_bindings_stack.top().end())
          {
            //mCRL2log(log::debug) << "  not found." << std::endl;
          }
          else
          {
            size_t index = f_it->second;
            //mCRL2log(log::debug) << "  index = " << index << " (" << print_index(index) << ")" << std::endl;
            //mCRL2log(log::debug) << "    => (" << print_index(index) << ") affects (" << print_index(offset) << ")" << std::endl;
            boost::add_edge(offset, index, G); // offset depends on index
          }
        }
      }
      offset++;
    }
    variable_bindings_stack.push(variable_bindings);
  }

  // two-pass is required for mes equations (all variables and parameters need to be known before
  // traversing the right hand sides.
  void register_mes_variable(const core::identifier_string& name, const data::variable_list& parameters)
  {
    propvars.push_back(name);
    propvar_offset_start[name] = offset; // start offset
    offset = propvar_offset_end[name] = offset + parameters.size(); // end offset
    for (auto p: parameters)
    {
      m_predicate_variables.push_back(p);
    }
  }

  void push_mes_variable(const core::identifier_string& name, const data::variable_list& parameters)
  {
    //mCRL2log(log::debug) << "Visit fixpoint formula: " << core::pp(name) << std::endl;
    variable_stack.push(name);
    size_t o = propvar_offset_start[name];
    std::map<data::variable, size_t> variable_bindings;
    if (!variable_bindings_stack.empty())
    {
      variable_bindings = variable_bindings_stack.top();
    }
    for (auto p: parameters)
    {
      variable_bindings[p] = o;
      o++;
    }
    variable_bindings_stack.push(variable_bindings);
  }

  void pop_variable()
  {
    variable_stack.pop();
    variable_bindings_stack.pop();
  }

  void enter(const state_formulas::nu& x)
  {
    push_variable(x.name(), x.assignments());
  }

  void leave(const state_formulas::nu& x)
  {
    pop_variable();
  }

  void enter(const state_formulas::mu& x)
  {
    push_variable(x.name(), x.assignments());
  }

  void leave(const state_formulas::mu& x)
  {
    pop_variable();
  }

  void enter(const state_formulas::must& x)
  {
    //mCRL2log(log::debug) << "Visit must: " << pp(x) << std::endl;
    regular_formulas::regular_formula x_formula = x.formula();
    std::set<data::variable> free_variables = regular_formulas::find_free_variables(x_formula);
    visit_free_variables(free_variables);
  }

  void enter(const state_formulas::may& x)
  {
    //mCRL2log(log::debug) << "Visit may: " << pp(x) << std::endl;
    regular_formulas::regular_formula x_formula = x.formula();
    std::set<data::variable> free_variables = regular_formulas::find_free_variables(x_formula);
    visit_free_variables(free_variables);
  }

  void apply(const modal_equation& x)
  {
    push_mes_variable(x.name(), x.parameters());
    (*this).apply(x.formula());
    pop_variable();
  }

  void apply(const modal_equation_system& x)
  {
    for(auto e: x.equations())
    {
      register_mes_variable(e.name(), e.parameters());
    }
    for(auto e: x.equations())
    {
      (*this).apply(e);
    }
  }

};

template <template <class> class Traverser>
struct apply_parelm_traverser: public Traverser<apply_parelm_traverser<Traverser> >
{
  typedef Traverser<apply_parelm_traverser<Traverser> > super;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::dependency_graph;
  using super::markings;
  using super::variables;
  using super::offset_start;
  using super::offset_end;

  apply_parelm_traverser(size_t N)
    : super(N)
  {}

#ifdef BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif
};


template <typename Derived>
struct count_variable_traverser: public state_formulas::state_formula_traverser<Derived>
{
  typedef state_formulas::state_formula_traverser<Derived> super;

  using super::enter;
  using super::leave;
  using super::apply;

#if BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif

  std::vector<core::identifier_string> m_variables;
  size_t m_parameter_count;

  count_variable_traverser() :
    m_parameter_count(0)
  {}

  const std::vector<core::identifier_string>& variables()
  {
    return m_variables;
  }

  size_t parameter_count()
  {
    return m_parameter_count;
  }

  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  void visit_variable(const core::identifier_string& name, const data::assignment_list& assignments)
  {
    m_variables.push_back(name);
    m_parameter_count += assignments.size();
  }

  void visit_variable(const core::identifier_string& name, const data::variable_list& parameters)
  {
    m_variables.push_back(name);
    m_parameter_count += parameters.size();
  }

  void enter(const state_formulas::nu& x)
  {
    visit_variable(x.name(), x.assignments());
  }

  void enter(const state_formulas::mu& x)
  {
    visit_variable(x.name(), x.assignments());
  }

  void apply(const modal_equation& x)
  {
    visit_variable(x.name(), x.parameters());
  }

  void apply(const modal_equation_system& x)
  {
    for(auto e: x.equations())
    {
      (*this).apply(e);
    }
  }
};

template <template <class> class Traverser>
struct apply_count_variable_traverser: public Traverser<apply_count_variable_traverser<Traverser> >
{
  typedef Traverser<apply_count_variable_traverser<Traverser> > super;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::variables;
  using super::parameter_count;

  apply_count_variable_traverser()
    : super()
  {}

#ifdef BOOST_MSVC
#include "mcrl2/core/detail/traverser_msvc.inc.h"
#endif
};


/// \brief Algorithm class for the parelm algorithm
class state_formula_parelm_algorithm
{
  protected:
    /// \brief The graph type of the dependency graph
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> graph;

    /// \brief The vertex type of the dependency graph
    typedef boost::graph_traits<graph>::vertex_descriptor vertex_descriptor;

    /// \brief The edge type of the dependency graph
    typedef boost::graph_traits<graph>::edge_descriptor edge_descriptor;

  public:

    /// \brief Runs the parelm algorithm. The state formula \p is modified by the algorithm
    /// \param x A state_formula
    template <typename T>
    void run(T& x) const
    {
      mCRL2log(log::debug) << "traversing formula: " << pp(x) << std::endl << std::endl;

      apply_count_variable_traverser<count_variable_traverser> vt;
      vt.apply(x);
      size_t N = vt.parameter_count();

      apply_parelm_traverser<parelm_traverser> f(N);
      f.apply(x);

      mCRL2log(log::debug) << "done traversing." << std::endl;

      graph G = f.dependency_graph();

      // the initial set v of significant variables
      std::set<size_t> v = f.markings();

      // compute the indices s of the parameters that need to be removed
      std::vector<size_t> r = utilities::reachable_nodes(G, v.begin(), v.end());
      std::sort(r.begin(), r.end());
      std::vector<size_t> q(N);
      utilities::detail::iota(q.begin(), q.end(), 0);
      std::vector<size_t> s;
      std::set_difference(q.begin(), q.end(), r.begin(), r.end(), std::back_inserter(s));

      // create a map that specifies the parameters that need to be removed
      std::map<core::identifier_string, std::vector<size_t> > removals;
      std::vector<size_t>::iterator sfirst = s.begin();
      for(auto it = f.variables().begin(); it != f.variables().end(); ++it)
      {
        core::identifier_string var = (*it);
        size_t offset_start = f.offset_start()[var];
        size_t offset_end = f.offset_end()[var];
        std::vector<size_t> offsets;
        while(sfirst != s.end() && (*sfirst) <= offset_end)
        {
          offsets.push_back((*sfirst) - offset_start);
          ++sfirst;
        }
        if (!offsets.empty())
        {
          removals[var] = offsets;
        }
      }

      // print debug output
      if (mCRL2logEnabled(log::debug))
      {
        mCRL2log(log::debug) << std::endl << "influential parameters:" << std::endl;
        for (std::set<size_t>::iterator i = v.begin(); i != v.end(); ++i)
        {
          //core::identifier_string X1 = find_predicate_variable(p, *i);
          data::variable v1 = f.predicate_variables()[*i];
          //mCRL2log(log::debug) << "(" + core::pp(X1) + ", " + data::pp(v1) + ")\n";
          mCRL2log(log::debug) << "(" << (*i) << ")" << std::endl;
          mCRL2log(log::debug) << "(" << pp(v1) << ")" << std::endl;
        }
        mCRL2log(log::debug) << std::endl << "dependencies:" << std::endl;
        typedef boost::graph_traits<graph>::edge_iterator edge_iterator;
        std::pair<edge_iterator, edge_iterator> e = edges(G);
        edge_iterator first = e.first;
        edge_iterator last  = e.second;
        for (; first != last; ++first)
        {
          edge_descriptor e = *first;
          size_t i1 = boost::source(e, G);
          size_t i2 = boost::target(e, G);
          mCRL2log(log::debug) << i1 << " -> " << i2 << std::endl;

        }
      }

      // print verbose output
      if (mCRL2logEnabled(log::verbose))
      {
        mCRL2log(log::verbose) << std::endl << "removing the following parameters:" << std::endl;
        for (std::map<core::identifier_string, std::vector<size_t> >::const_iterator i = removals.begin(); i != removals.end(); ++i)
        {
          core::identifier_string X1 = i->first;

          for (std::vector<size_t>::const_iterator j = (i->second).begin(); j != (i->second).end(); ++j)
          {
            data::variable v1 = f.predicate_variables()[*j + f.offset_start()[X1]];
            mCRL2log(log::verbose) << "(" << core::pp(X1) << ", " << (*j) << " | " << pp(v1) << ": " << pp(v1.sort()) << ")" << std::endl;
          }
        }
      }

      state_formulas::algorithms::remove_parameters(x, removals);
    }
};

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_PARELM_H
