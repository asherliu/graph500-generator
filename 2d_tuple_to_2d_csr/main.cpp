#include "translator.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>


int main(int argc, char* argv[]) {


	index_t *ranger_beg_pos, *ranger_adj_card;
	data_t *ranger_adj_list;

	translator<data_t, index_t>
	(argc, argv, ranger_adj_list, ranger_beg_pos, 
	ranger_adj_card);

	return 0;
}
