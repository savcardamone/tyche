/**
 * @file multiprocess_communications.cpp
 * @author Salvatore Cardamone
 * @brief Communications construction and management.
 */
#include <set>
#include <vector>
#include "multiprocess/multiprocess_communications.hpp"

namespace tycheplusplus {

/**
 * @brief Class constructor. Initialises communicator as world.
 */
MultiProcessCommunications::MultiProcessCommunications()
    : group_id_(0) {}

/**
 * @brief Class constructor. Explicitly initialise communicator with argument.
 * @param comm Communicator to initialise with
 */
MultiProcessCommunications::MultiProcessCommunications(
    boost::mpi::communicator& comm) : group_id_(0), comm_(comm) {}

/**
 * @brief Partition the world communicator into groups based on input index.
 *        We're also able to initialise the group leader communicator with this.
 * @param group_id Index of the group to which the process belongs.
 */
void MultiProcessCommunications::Partition(const unsigned int& group_id) {

  // Split communicator into groups based on input index
  group_id_ = group_id;
  group_comm_ = boost::mpi::communicator(comm_.split(group_id));

}

/**
 * @brief Partition the world communicator into groups based on nodes upon which
*         processes reside.
 */
void MultiProcessCommunications::PartitionNodes() {

  // Node ID the process resides on
  unsigned int node_id = 0;

  // If this is the master process, gather the processor name from all other
  // processes
  if (comm_.rank() == 0) {

    // Gather all process processor names
    std::vector<std::string> processor_names;
    boost::mpi::gather(
        comm_, boost::mpi::environment::processor_name(), processor_names, 0);

    // Initialise the set with the vector. Removes duplicates
    std::set<std::string> set_names(
        processor_names.begin(), processor_names.end());

    // Get the node index of the master process. Note that this isn't
    // necessarily zero
    node_id = std::distance(
        set_names.begin(), set_names.find(processor_names[0]));

    // Loop through the processor name for each process and locate its "index"
    // from the unique set
    for (int i_proc = 1; i_proc < comm_.size(); ++i_proc) {
      int message = std::distance(
          set_names.begin(), set_names.find(processor_names[i_proc]));
      comm_.send(i_proc, 0, message);
    }

  // If this is a slave process, send the process name to the master slave
  } else {

    // Send processor name to master process
    boost::mpi::gather(comm_, boost::mpi::environment::processor_name(), 0);
    // Receive a node index from the master process
    comm_.recv(0, 0, node_id);

  }

  // Perform the communicator partitioning
  Partition(node_id);

}

/**
 * @brief Print some generic information about the communicators in the object.
 * @param stream Output stream to dump information to.
 */
void MultiProcessCommunications::Print(std::ostream& stream) const {

  stream << " === MultiProcessCommunications" << std::endl ;
  stream << "     " << comm_.rank() << " of " << comm_.size()
         << " processes in world." << std::endl ;
  stream << "     " << group_comm_.rank() << " of " << group_comm_.size()
         << " processes in group." << std::endl ;

}

} /* namespace tycheplusplus */
