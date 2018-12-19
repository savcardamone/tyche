/**
 * @file output_manager.hpp
 * @author Salvatore Cardamone
 * @brief Output manager class for streams.
 */
#ifndef __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP
#define __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP

#include <memory>
#include <fstream>
#include <iostream>
#include "multiprocess/multiprocess_base.hpp"

namespace tycheplusplus {

/**
 * @class OutputManager
 * @brief Manage output streams. Primarily encapsulated so we can control
 *        processes we don't want to print to stdout or files.
 */
class OutputManager : public MultiProcessBase {

public:
  OutputManager(boost::mpi::communicator& comm, std::streambuf* buffer);
  std::ostream& GetStream() { return outstream_; };

private:
  void SwitchBuffers();
  std::ostream outstream_, nullstream_ ;

};

/**
 * @brief << Operator overload for the OutputManager. Streams a value into the
 *        std::ostream within.
 * @param manager OutputManager object.
 * @param val Value to place on the std::ostream.
 * @retval The OutputManager object. Effectively turns a multi-valued stream
 *         into a recursion.
 */
template<class T>
inline OutputManager& operator<<(OutputManager& manager, const T& val) {

  manager.GetStream() << val ;
  return manager ;

}

/**
 * @brief << overload for the std::endl, which is a template function rather
 *        than a value. Won't default to above overload with std::endl, will
 *        take something from the standard library instead.
 * @ref https://stackoverflow.com/a/38077832
 */
inline OutputManager& operator<<(
    OutputManager& manager, std::ostream& (*val)(std::ostream &)) {

  manager.GetStream() << val ;
  return manager ;

}

} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_OUTPUT_MANAGER_HPP */
