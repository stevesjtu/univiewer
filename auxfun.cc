#include "auxfunc.h"

void argParser(const int& argc, char* argv[], 
                        std::vector<std::string> &simple_output_result,
                        std::vector<std::string> &modelFiles,
											  std::vector<std::string> &dispFiles,
											  std::vector<std::string> &contFiles,
											  std::vector<std::string> &nodeDatafiles)
{
	if (argc == 1) {
		std::cout << "Type 'Univiewer /h' for more help." << std::endl;
		exit(0);
	}
	std::vector<std::string> *ptr_vector_name = NULL;
	for (int i = 1; i<argc; ++i) {

#ifdef __APPLE__
		if (*argv[i] == '-') {
#else
		if ((*argv[i] == '-') || (*argv[i] == '/')) {
#endif
			switch (*(argv[i] + 1)) {
			case 'm':
				ptr_vector_name = &modelFiles;
				break;
			case 'o':
				ptr_vector_name = &dispFiles;
				break;
			case 'r':
				ptr_vector_name = &simple_output_result;
				break;
			case 's':
				ptr_vector_name = &nodeDatafiles;
				break;
			case 'c':
				ptr_vector_name = &contFiles;
				break;
			case 'h':
				std::cout << std::endl;
				std::cout << "Usage: Univiewer /m file1.xml file2.xml ... fileN.xml /o disp.dat" << std::endl;
				std::cout << std::endl;
				std::cout << "The 'fileN.xml' is a model file that is compatible with VTK API, and other programs. It is a text file using *.vtu format, but you can also add your own data." << std::endl;
				std::cout << "The 'disp.dat' is a data file that contains the displacements of nodes for all the models(file1.xml, file2.xml and so on). It is a binary file with all the data as the type of double, the detail data sequence is as follows:" << std::endl << std::endl;
				std::cout << "#####################################################################################################" << std::endl;
				std::cout << "time1 \n node_1_disp_x node_1_disp_y node_1_disp_z \n node_2_disp_x node_2_disp_y node_2_disp_z \n ... \n node_n_disp_x node_n_disp_y node_n_disp_z" << std::endl;
				std::cout << "time2 \n ......" << std::endl;
				std::cout << "#####################################################################################################" << std::endl << std::endl;
				std::cout << "NOTE: Binary files are not same as txt files which contains symbols like \\n \\t, and all the data saved as String. A Binary file contains data saved as its own type without any symbols like \\n \\t." << std::endl;
				exit(0);
				break;
			default:
				break;
			}
			continue;
		}
		ptr_vector_name->push_back(argv[i]);
		}
}
