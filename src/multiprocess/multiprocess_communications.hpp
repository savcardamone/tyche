/**
 * @file multiprocess_communications.hpp
 * @author Salvatore Cardamone
 * @brief Multiprocess class which allows us to construct and manage
 *        communicators.
 */
#ifndef __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP
#define __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP

#include <ostream>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace tycheplusplus {
  
/**
 * @class MultiProcessCommunications
 * @brief Multiprocess communications construction and management.
 */
class MultiProcessCommunications {

public:
  MultiProcessCommunications();

  const boost::mpi::communicator& WorldComm() const { return world_comm_; }
  boost::mpi::communicator& WorldComm() { return world_comm_; }
  
  const boost::mpi::communicator& NodeComm() const { return node_comm_; }
  boost::mpi::communicator& NodeComm() { return node_comm_; }

  inline int WorldSize() const { return WorldComm().size(); }
  inline int NodeSize() const { return NodeComm().size(); }

  inline int WorldRank() const { return WorldComm().rank(); }
  inline int NodeRank() const { return NodeComm().rank(); }

  inline int NodeID() const { return node_id_; }

  inline bool IsWorldMaster() const { return WorldRank() == 0 ? true : false; }
  inline bool IsNodeMaster() const { return NodeRank() == 0 ? true : false; }

  friend std::ostream& operator<<(std::ostream& os,
				  const MultiProcessCommunications& m);
  
private:
  int node_id_;
  boost::mpi::environment env_;
  boost::mpi::communicator world_comm_, node_comm_;

} ;
  
} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_MULTIPROCESS_COMMUNICATIONS_HPP */
