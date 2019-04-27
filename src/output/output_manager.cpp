/**
 * @file output_manager.cpp
 * @author Salvatore Cardamone
 * @brief Output manager implementations.
 */
#include <iostream>
#include "output/output_manager.hpp"

namespace tycheplusplus {

/**
 * @brief Class constructor. Initialise the output streams and set the active stream
 *        as requested. Note that the active stream of a process which isn't the master 
 *        of the communicator defaults to the nullstream. 
 * @param comm Reference to MPI communicator.
 * @param name Unique name for the OutputManager to initialise a file with when/if required.
 * @param stream Type of the stream the communicator master initialises its active stream with.
 */
OutputManager::OutputManager(boost::mpi::communicator& comm, std::string& name, StreamTypes stream)
  : MultiProcessBase(comm), name_(name), term_stream_(std::cout), null_stream_(nullptr), active_stream_(nullptr) {

  // If process isn't communicator master, redirect its output to
  // the nullstream
  if (IsMaster() == true) {
    SwitchStream(stream);
  } else {
    SwitchStream(StreamTypes::Null);
  }

}

/**
 * @brief Class destructor. Check whether we need to close the object's file
 *        member if it was opened.
 */
OutputManager::~OutputManager() {

  // Close the file stream if it's open to flush the buffer
  if (file_stream_.is_open()) file_stream_.close();
    
}
  
/**
 * @brief Switch the active stream to the requested stream.
 * @param stream The stream type to switch the active stream to.
 */
void OutputManager::SwitchStream(const StreamTypes stream) {

  switch(stream) {

    // Output is directed to terminal
    case StreamTypes::Terminal:
      active_stream_.rdbuf(term_stream_.rdbuf());
      break;

    // Output is directed to file
    case StreamTypes::File:
      // If we haven't opened the file stream yet, then do so now
      if (file_stream_.is_open() == false) {
        std::string filename(name_ + ".log");
        file_stream_.open(filename);
      }
      active_stream_.rdbuf(file_stream_.rdbuf());
      break;

    // Output is directed to somewhere we don't care about
    case StreamTypes::Null:
      active_stream_.rdbuf(null_stream_.rdbuf());
      break;

    default:
      throw std::invalid_argument("StreamType is not supported.");
      break;

  }

}
  
} /* namespace tycheplusplus */
