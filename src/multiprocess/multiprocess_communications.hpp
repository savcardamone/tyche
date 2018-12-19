/**
 * @file multiprocess_communications.hpp
 * @author Salvatore Cardamone
 * @brief Multiprocess class which allows us to construct and manage
 *        communicators.
 */
#ifndef __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP
#define __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP

#include <iostream>
#include <fstream>
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>

namespace tycheplusplus {

/**
 * @class MultiProcessCommunications
 * @brief Multiprocess communications construction and management.
 */
class MultiProcessCommunications {

public:
  MultiProcessCommunications();
  MultiProcessCommunications(boost::mpi::communicator& comm);

  void PartitionNodes();
  void Print(std::ostream& stream) const;

  boost::mpi::communicator& GetCommunicator() {
    return comm_ ;
  }

private:
  unsigned int group_id_;
  boost::mpi::environment env_;
  boost::mpi::communicator comm_, group_comm_;

  void Partition(const unsigned int& group_id);

} ;

} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP */
