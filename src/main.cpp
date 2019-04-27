/**
 * @file main.cpp
 * @author Salvatore Cardamone
 * @brief Entry point for tyche++.
 */
#include "multiprocess/multiprocess_communications.hpp"
#include "output/output_manager.hpp"

/**
 * @brief Entry routine for tyche++.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @retval EXIT_SUCCESS on success, else EXIT_FAILURE.
 */
int main(int argc, char* argv[]) {

  tycheplusplus::MultiProcessCommunications comms;

  std::string world_output_name(
     "World_Proc" + std::to_string(comms.WorldRank())
  );
  tycheplusplus::OutputManager world_output(
  	comms.WorldComm(), world_output_name, tycheplusplus::OutputManager::StreamTypes::Terminal
  );
  world_output << comms;

  std::string node_output_name(
     "Node" + std::to_string(comms.NodeID()) + "_Proc" + std::to_string(comms.WorldRank())
  );
  tycheplusplus::OutputManager node_output(
  	comms.NodeComm(), node_output_name, tycheplusplus::OutputManager::StreamTypes::File
  );
  node_output << comms;

}
