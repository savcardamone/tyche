/**
 * @file output_manager.cpp
 * @author Salvatore Cardamone
 * @brief Output manager implementations.
 */
#include "output/output_manager.hpp"

namespace tycheplusplus {

/**
 * @brief Class constructor. Initialise the output stream with argument and the
 *        temp stream with the nullstream.
 * @param comm Shared pointer to MPI communicator.
 * @param stream Shared pointer to output stream.
 */
OutputManager::OutputManager(
    boost::mpi::communicator& comm, std::streambuf* buffer)
    : MultiProcessBase(comm), outstream_(buffer), nullstream_(nullptr) {

  if (Rank() != 0) SwitchBuffers();

}

/**
 * @brief Switch the nullstream and output stream buffers. Either activates or
 *        deactivates the output stream, depending on the start state. If
 *        OutputManager was constructed with a nullptr, the this method does
 *        nothing.
 */
void OutputManager::SwitchBuffers() {

  // Temporary store the streambuf in a temporary
  auto temp_buffer = nullstream_.rdbuf();
  // Swap the nullstream and outstream streambufs
  nullstream_.rdbuf(outstream_.rdbuf()); outstream_.rdbuf(temp_buffer);

}

} /* namespace tycheplusplus */
