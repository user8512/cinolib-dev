#include "patch.h"

int main() {
#ifdef TEST
	cinolib::Patch patch(root + "/input/block.mesh");
	std::cout << "已读取" << patch.getMesh()->vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(3);
#else
	cinolib::Patch patch(root + "/input/bunny_hex.mesh");
	std::cout << "已读取" << patch.getMesh()->vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(1);
#endif
}