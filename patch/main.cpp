#include "patch.h"


int main() {
#ifdef USE_CUDA
	if (!cuda::init()) {
		std::cout << "cuda init error" << std::endl;
	}
#endif
#ifdef TEST
	cinolib::Patch patch(root + "/input/block.mesh");
	std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(5);
#else
	cinolib::Patch patch(root + "/input/bunny_hex.mesh");
	std::cout << "已读取" << patch.getMesh().vector_polys().size() << "单元体网格" << std::endl;
	patch.subdiv(1);
#endif
}