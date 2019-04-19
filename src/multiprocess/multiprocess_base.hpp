/**
 * @file multiprocess_base.hpp
 * @author Salvatore Cardamone
 * @brief Base multiprocess class from which multiprocessing-aware classes
 *        can derive.
 */
#ifndef __TYCHEPLUSPLUS_MULTIPROCESS_BASE_HPP
#define __TYCHEPLUSPLUS_MULTIPROCESS_BASE_HPP

#include <boost/mpi/communicator.hpp>

namespace tycheplusplus {

/**
 * @class MultiProcessBase
 * @brief Base MPI class from which we can derive.
 */
class MultiProcessBase {

public:
  /**
   * @brief Class constructor. Initialise communicator.
   * @param comm_ Communicator for initialisation.
   */
  MultiProcessBase(boost::mpi::communicator& comm) : comm_(comm) {}

protected:
  /**
   * @brief Return the size of the communicator.
   * @retval Number of processes in the communicator.
   */
  inline int Size() const { return comm_.size(); }

  /**
   * @brief Return the rank of the process on the communicator.
   * @retval Rank of the process.
   */
  inline int Rank() const { return comm_.rank(); }

  /**
   * @brief Query whether the process in the communicator master.
   * @retval True if it's master (i.e. rank zero), false otherwise.
   */
  inline bool IsMaster() const { return comm_.rank() == 0 ? true : false; }

private:
  // Communicator the deriving object can use
  boost::mpi::communicator comm_;

};

} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_MULTIPROCESS_BASE_HPP */
