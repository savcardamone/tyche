/**
 * @file main.cpp
 * @author Salvatore Cardamone
 * @brief Entry point for tyche++.
 */
#include "multiprocess/multiprocess_communications.hpp"
#include "output/output_manager.hpp"
#include "input/input_manager.hpp"
#include "utilities/fixed_point.hpp"

/**
 * @brief Entry routine for tyche++.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @retval EXIT_SUCCESS on success, else EXIT_FAILURE.
 */
int main(int argc, char* argv[]) {

  tycheplusplus::MultiProcessCommunications comms;
  tycheplusplus::OutputManager output(comms.GetCommunicator(), std::cout.rdbuf());
  output << "tyche++ begins" << std::endl;

  comms.PartitionNodes();
  comms.Print(output.GetStream());
  
  tycheplusplus::InputManager input();

  output << "tyche++ ends." << std::endl;

}
