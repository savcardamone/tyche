/**
 * @file output_manager.hpp
 * @author Salvatore Cardamone
 * @brief Multiprocessing-aware output manager class.
 */
#ifndef __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP
#define __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP

#include <fstream>
#include "multiprocess/multiprocess_base.hpp"

namespace tycheplusplus {

/**
 * @class OutputManager
 * @brief Manage output streams. Encapsulates the streams we can output to and
 *        direct output to the appropriate stream.
 */
class OutputManager : public MultiProcessBase {

public:

  // Stream types we support with the OutputManager
  enum StreamTypes { 
    Terminal, File, Null
  };

  // Constructors and Destructors
  OutputManager(boost::mpi::communicator& comm, std::string& name, StreamTypes stream);
  ~OutputManager();
  
  /**
   * @brief Return a reference to the active stream.
   * @retval Reference to the object's active stream.
   */
  inline std::ostream& Stream() { return active_stream_; }

  void SwitchStream(const StreamTypes Stream);
  
  /**
   * @brief << operator overload for the OutputManager. Streams a value into the
   *        active stream. Note that we serialise over processes in case
   *        processes share an ostream (i.e. the terminal).
   * @param manager OutputManager reference.
   * @param val Value to place on the active stream.
   * @retval The OutputManager reference.
   */
  template<class T>
  friend OutputManager& operator<<(OutputManager& manager, const T& val) {

    // Go through processes one by one and print to the ostream one at a time
    for (int i_proc=0; i_proc<manager.Size(); ++i_proc) {
      if (i_proc == manager.Rank()) manager.Stream() << val;
      manager.Comm().barrier();
    }

    return manager ;

  }
  
private:
  // Name of the OutputManager. For multiprocessing, we can then get a unique filename
  // we can write to if output is directed to the file stream
  std::string name_;
  // Note that since std::ostream is not copy-constructible, we can only have a reference
  // to the terminal stream since it's initialised with std::cout
  std::ostream& term_stream_; 
  std::ostream null_stream_, active_stream_;
  std::ofstream file_stream_;

};

} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP */
