// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/modal_equation_system.h
/// \brief Add your file description here.

#ifndef MCRL2_MODAL_EQUATION_SYSTEM_H_
#define MCRL2_MODAL_EQUATION_SYSTEM_H_

#include <cassert>
#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/atermpp/aterm_list.h"
#include "mcrl2/core/detail/default_values.h"
#include "mcrl2/core/detail/soundness_checks.h"
#include "mcrl2/core/detail/function_symbols.h"
#include "mcrl2/data/data_specification.h"
#include "mcrl2/modal_formula/state_formula.h"

namespace mcrl2
{

namespace state_formulas
{

/// \brief A fixpoint symbol
class fixpoint_symbol: public atermpp::aterm_appl
{
  public:
    /// \brief Default constructor.
    fixpoint_symbol()
      : atermpp::aterm_appl(core::detail::default_values::FixPoint)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit fixpoint_symbol(const atermpp::aterm& term)
      : atermpp::aterm_appl(term)
    {
      assert(core::detail::check_rule_FixPoint(*this));
    }

    /// \brief Returns the mu symbol.
    /// \return The mu symbol.
    static fixpoint_symbol mu()
    {
      return fixpoint_symbol(atermpp::aterm_appl(core::detail::function_symbol_Mu()));
    }

    /// \brief Returns the nu symbol.
    /// \return The nu symbol.
    static fixpoint_symbol nu()
    {
      return fixpoint_symbol(atermpp::aterm_appl(core::detail::function_symbol_Nu()));
    }

    /// \brief Returns true if the symbol is mu.
    /// \return True if the symbol is mu.
    bool is_mu() const
    {
      return function() == core::detail::function_symbols::Mu;
    }

    /// \brief Returns true if the symbol is nu.
    /// \return True if the symbol is nu.
    bool is_nu() const
    {
      return function() == core::detail::function_symbols::Nu;
    }
};

// prototype declaration
std::string pp(const fixpoint_symbol& x);

/// \brief Outputs the object to a stream
/// \param out An output stream
/// \return The output stream
inline
std::ostream& operator<<(std::ostream& out, const fixpoint_symbol& x)
{
  return out << state_formulas::pp(x);
}



/// \brief A modal equation
class modal_equation
{
  protected:
    /// \brief The fixpoint symbol of the equation
    fixpoint_symbol m_symbol;

    /// \brief The variable name on the left hand side of the equation
    core::identifier_string m_name;

    /// \brief The data parameters of the variable
    data::variable_list m_parameters;

    /// \brief The expression on the right hand side of the equation
    state_formula m_formula;

  public:
    modal_equation()
    {}

    /// \brief Default constructor.
    modal_equation(const fixpoint_symbol& symbol, const core::identifier_string& name,
        const data::variable_list& parameters, const state_formula& formula) :
      m_symbol(symbol),
      m_name(name),
      m_parameters(parameters),
      m_formula(formula)
    {}

    /// \brief Copy constructor.
    modal_equation(const modal_equation& other) :
      m_symbol(other.m_symbol),
      m_name(other.m_name),
      m_parameters(other.m_parameters),
      m_formula(other.m_formula)
    {}

    fixpoint_symbol& symbol()
    {
      return m_symbol;
    }

    const fixpoint_symbol& symbol() const
    {
      return m_symbol;
    }

    core::identifier_string& name()
    {
      return m_name;
    }

    const core::identifier_string& name() const
    {
      return m_name;
    }

    data::variable_list& parameters()
    {
      return m_parameters;
    }

    const data::variable_list& parameters() const
    {
      return m_parameters;
    }

    void set_parameters(const data::variable_list& parameters)
    {
      m_parameters = parameters;
    }

    state_formula& formula()
    {
      return m_formula;
    }

    const state_formula& formula() const
    {
      return m_formula;
    }

    void set_formula(const state_formula& formula)
    {
      m_formula = formula;
    }
};

// prototype declaration
std::string pp(const modal_equation& x);

/// \brief Outputs the object to a stream
/// \param out An output stream
/// \return The output stream
inline
std::ostream& operator<<(std::ostream& out, const modal_equation& x)
{
  return out << state_formulas::pp(x);
}



/// \brief A modal equation system
class modal_equation_system
{
  protected:
    variable m_initial_state;
    std::vector<modal_equation> m_equations;
    std::map<core::identifier_string, modal_equation> m_equation_map;

  public:
    /// \brief Default constructor.
    modal_equation_system()
    {}

    /// \brief Copy constructor.
    modal_equation_system(const modal_equation_system& other) :
      m_initial_state(other.m_initial_state),
      m_equations(other.m_equations),
      m_equation_map(other.m_equation_map)
    {}

    modal_equation_system(const variable& initial_state, const std::vector<modal_equation>& equations) :
      m_initial_state(initial_state),
      m_equations(equations)
    {
      for(auto e : m_equations)
      {
        if (m_equation_map.find(e.name()) != m_equation_map.end())
        {
          throw std::runtime_error("Equation names are not unique.");
        }
        m_equation_map[e.name()] = e;
      }
    }

    void set_initial_state(const variable& initial_state)
    {
      m_initial_state = initial_state;
    }

    variable& initial_state()
    {
      return m_initial_state;
    }

    const variable& initial_state() const
    {
      return m_initial_state;
    }

    std::vector<modal_equation>& equations()
    {
      return m_equations;
    }

    const std::vector<modal_equation>& equations() const
    {
      return m_equations;
    }

    std::map<core::identifier_string, modal_equation>& equation_map()
    {
      return m_equation_map;
    }

    const std::map<core::identifier_string, modal_equation>& equation_map() const
    {
      return m_equation_map;
    }

    void set_equations(const std::vector<modal_equation>& equations)
    {
      m_equations = equations;
      m_equation_map.clear();
      for(auto e : m_equations)
      {
        if (m_equation_map.find(e.name()) != m_equation_map.end())
        {
          throw std::runtime_error("Equation names are not unique.");
        }
        m_equation_map[e.name()] = e;
      }
    }

};

// prototype declaration
std::string pp(const modal_equation_system& x);

/// \brief Outputs the object to a stream
/// \param out An output stream
/// \return The output stream
inline
std::ostream& operator<<(std::ostream& out, const modal_equation_system& x)
{
  return out << state_formulas::pp(x);
}

} // namespace state_formulas

} // namespace mcrl2

#endif /* MCRL2_MODAL_EQUATION_SYSTEM_H_ */
