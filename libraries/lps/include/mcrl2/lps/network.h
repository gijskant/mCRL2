// Author(s): Gijs Kant
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/network.h
/// \brief The class network.

#ifndef MCRL2_LPS_NETWORK_H_
#define MCRL2_LPS_NETWORK_H_

#include "mcrl2/lps/specification.h"
#include "mcrl2/lps/synchronization_vector.h"

namespace mcrl2
{

namespace lps
{

class network; // prototype declaration

///////////////////////////////////////////////////////////////////////////////
// synchronization_vector
/// \brief synchronization vector.
class network
{
  protected:
    /// \brief The synchronization vector stored as a map from the resulting action label
    /// to a vector of vectors that contain one action label for every LPS in the system.
    std::vector<std::string> m_lps_filenames;
    lps::synchronization_vector m_synchronization_vector;
    process::action_label_list m_action_labels;

  public:
    /// \brief Constructor.
    network();

    /// \brief Copy constructor.
    network(const network &other);

    /// \brief Constructor.
    network(std::vector<std::string> lps_filenames,
            lps::synchronization_vector synchronization_vector,
            process::action_label_list action_labels);

    /// \brief Reads the network from a stream.
    /// \param in A stream
    /// Processes are represented by separate files as specified in the file filename.
    void read(std::istream& in);

    /// \brief Reads the network from file.
    /// \param filename A string
    /// If filename is nonempty, input is read from the file named filename.
    /// If filename is empty, input is read from standard input.
    /// Processes are represented by separate files as specified in the file filename.
    void load(const std::string& filename);

    /// \brief Writes the network to a stream.
    /// \param out A stream
    /// The processes are not written, only the filenames of the original processes.
    void write(std::ostream& out) const;

    /// \brief Writes the network to file.
    /// \param filename A string
    /// If filename is nonempty, output is written to the file named filename.
    /// If filename is empty, output is written to stdout.
    /// The processes are not written to file, only the filenames of the original processes.
    void save(const std::string& filename) const;

    /// \brief Destructor
    ~network();

    /// \brief Returns the vector of LPS filenames.
    /// \return The vector of filenames.
    const std::vector<std::string>& lps_filenames() const;

    /// \brief Returns the vector of LPS filenames.
    /// \return The vector of filenames.
    std::vector<std::string>& lps_filenames();

    /// \brief Returns the synchronization vector.
    /// \return The synchronization vector.
    const lps::synchronization_vector& synchronization_vector() const;

    /// \brief Returns the synchronization vector.
    /// \return The synchronization vector.
    lps::synchronization_vector& synchronization_vector();

    /// \brief Returns the actions labels.
    /// \return The action label list.
    const process::action_label_list& action_labels() const;

    /// \brief Returns the actions labels.
    /// \return The action label list.
    process::action_label_list& action_labels();

};

// template function overloads
std::string pp(const network& x);

} // namespace lps

} // namespace mcrl2


#endif /* MCRL2_LPS_NETWORK_H_ */
