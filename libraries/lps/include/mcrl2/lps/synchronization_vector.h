// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/synchronization_vector.h
/// \brief The class synchronization_vector.

#ifndef MCRL2_LPS_SYNCHRONIZATION_VECTOR_H
#define MCRL2_LPS_SYNCHRONIZATION_VECTOR_H

#include "mcrl2/process/action_label.h"

namespace mcrl2
{

namespace lps
{

typedef std::pair<std::vector<std::string>, process::action_label > synchronization_vector_element;
typedef std::vector<synchronization_vector_element> synchronization_vector_type;

static const std::string inactive_label = "inactive";

class synchronization_vector; // prototype declaration

///////////////////////////////////////////////////////////////////////////////
// synchronization_vector
/// \brief synchronization vector.
class synchronization_vector
{

  protected:
    /// \brief The synchronization vector stored as a map from the resulting action label
    /// to a vector of vectors that contain one action label for every LPS in the system.
    synchronization_vector_type m_synchronization_vector;
    /// \brief the set of action labels in the synchronization vector.
    std::set<process::action_label> m_action_labels;
    /// \brief the length of the vectors, i.e., the number of components.
    size_t m_length;

  public:
    /// \brief Constructor.
    synchronization_vector();

    /// \brief Copy constructor.
    synchronization_vector(const synchronization_vector &other);

    /// \brief Constructor.
    synchronization_vector(int length);

    /// \brief Constructor.
    synchronization_vector(synchronization_vector_type synchronization_vector,
        std::set<process::action_label> action_labels,
        int length);

    void normalize(const data::data_specification& dataspec);

    /// \brief Reads the synchronization vector from stream.
    /// \param in An input stream
    void read(std::istream& in);

    /// \brief Reads the synchronization vector from file.
    /// \param filename A string
    /// If filename is nonempty, input is read from the file named filename.
    /// If filename is empty, input is read from standard input.
    void load(const std::string& filename);

    /// \brief Writes the synchronization vector to a stream.
    /// \param out An output stream
    void write(std::ostream& out) const;

    /// \brief Writes the synchronization vector to file.
    /// \param filename A string
    /// If filename is nonempty, output is written to the file named filename.
    /// If filename is empty, output is written to stdout.
    void save(const std::string& filename) const;

    /// \brief Destructor
    ~synchronization_vector();

    void add_vector(std::vector<std::string> actions, process::action_label a);

    const synchronization_vector_type& vector() const;

    const std::set<process::action_label>& action_labels() const;

    size_t length() const;

};

// template function overloads
std::string pp(const synchronization_vector& x);
std::string pp(const synchronization_vector_element& x);

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_SYNCHRONIZATION_VECTOR_H
