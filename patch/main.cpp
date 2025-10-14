#include "patch.h"

std::ofstream outlog(root + "/output/test/log.txt");

int main() {
#ifdef TEST
	cinolib::Patch patch(root + "/input/rockerarm_subdiv_3.mesh");
	//std::cout << "已读取" << patch.getMesh()->vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(1);
#else
	cinolib::Patch patch(root + "/input/rockerarm_subdiv_3.mesh");
	//std::cout << "已读取" << patch.getMesh()->vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(1);
#endif
}