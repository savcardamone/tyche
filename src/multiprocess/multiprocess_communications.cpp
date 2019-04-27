/**
 * @file multiprocess_communications.cpp
 * @author Salvatore Cardamone
 * @brief Communications construction and management.
 */
#include <set>
#include <vector>
#include <boost/mpi/collectives.hpp>

#include "multiprocess/multiprocess_communications.hpp"

namespace tycheplusplus {

/**
 * @brief Class constructor. Initialises communicator as world and constructs
 *        a separate communicator for all processes on a node.
 */
MultiProcessCommunications::MultiProcessCommunications() {

  static const int master_proc_id = 0;
  // Node ID the process resides on
  node_id_ = 0;

  // If this is the master process, gather the processor name from all other
  // processes
  if (world_comm_.rank() == 0) {

    // Gather all process processor names
    std::vector<std::string> processor_names;
    boost::mpi::gather(
        world_comm_, boost::mpi::environment::processor_name(),
	      processor_names, master_proc_id
    );

    // Initialise the set with the vector. Removes duplicates
    std::set<std::string> set_names(
        processor_names.begin(), processor_names.end());

    // Get the node index of the master process. Note that this isn't
    // necessarily zero
    node_id_ = std::distance(
        set_names.begin(), set_names.find(processor_names[master_proc_id]));

    // Loop through the processor name for each process and locate its "index"
    // from the unique set
    for (int i_proc = master_proc_id+1; i_proc < world_comm_.size(); ++i_proc) {
      int message = std::distance(
          set_names.begin(), set_names.find(processor_names[i_proc]));
      world_comm_.send(i_proc, 0, message);
    }

  // If this is a slave process, send the process name to the master slave
  } else {

    // Send processor name to master process
    boost::mpi::gather(
        world_comm_, boost::mpi::environment::processor_name(), 0
    );

    // Receive a node index from the master process
    world_comm_.recv(master_proc_id, 0, node_id_);

  }

  // Split communicator into groups based on input index
  node_comm_ = boost::mpi::communicator(world_comm_.split(node_id_));
  
}

/**
 * @brief Dump some information about the multiprocess communications to an
 *        ostream.
 * @param os An output stream iterator.
 * @param m MultiProcessCommunications object we're dumping.
 * @retval An output stream iterator.
 */
std::ostream& operator<<(std::ostream& os, const MultiProcessCommunications& m) {

  os << m.WorldRank() << " of " << m.WorldSize() << " processes in world." << std::endl ;
  os << m.NodeRank()  << " of " << m.NodeSize()  << " processes on node."  << std::endl ;

  return os;
  
}

} /* namespace tycheplusplus */
