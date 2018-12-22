/**
 * @file input_manager.hpp
 * @author Salvatore Cardamone
 * @brief Base input manager class.
 */
#ifndef __TYCHEPLUSPLUS_INPUT_MANAGER_HPP
#define __TYCHEPLUSPLUS_INPUT_MANAGER_HPP

#include <libxml++/libxml++.h>

namespace tycheplusplus {

/**
 * @class InputManager
 * @brief
 */
class InputManager {

public:
  InputManager() {
    // Set the global C++ locale to the user-specified locale. Then we can
    // hopefully use std::cout with UTF-8, via Glib::ustring, without exceptions.
    std::locale::global(std::locale(""));
    // xml documents will be validated by default
    parser_.set_validate();
  }

  void ParseFile(std::string filename) {
    parser_.parse_file(filename);
  }
  
private:
  xmlpp::DomParser parser_;
  
};

} /* namespace tycheplusplus */

#endif /* #ifndef __TYCHEPLUSPLUS_INTPUT_MANAGER_HPP */
