#include "SimFramework.h"
#include <QtWidgets/QApplication>
#include "Commom/PolynomialSolver.h"

bool detectCudaDevice(double& gpu_mem);

int main(int argc, char *argv[])
{
	double gpu_mem = 0;
	if (!detectCudaDevice(gpu_mem))
		return 0;

	QApplication a(argc, argv);
	SimFramework w;
	BaseMainWidget* simWidget = new BaseMainWidget();
	simWidget->connectUI(&w);

	BaseSimulator* sim = new BaseSimulator();
	simWidget->getCenteralScene()->bindSimulator(sim);

	w.setWindowTitle("Sim Framework");
	w.setWindowIcon(QIcon("./window_ico.ico"));
	w.showMaximized();
	return a.exec();
}

bool detectCudaDevice(double& gpu_mem)
{
	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
	if (error_id != cudaSuccess) {
		printf("cudaGetDeviceCount returned %d\n-> %s\n",
			static_cast<int>(error_id), cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		system("pause");
		return false;
	}

	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount == 0) {
		printf("There are no available device(s) that support CUDA\n");
		system("pause");
		return false;
	}
	else {

		int dev, driverVersion = 0, runtimeVersion = 0;

		for (dev = 0; dev < deviceCount; ++dev) {
			cudaSetDevice(dev);
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);
			printf("Detected a CUDA Capable device(s): \"%s\" \n", deviceProp.name);
			// Console log
			cudaDriverGetVersion(&driverVersion);
			cudaRuntimeGetVersion(&runtimeVersion);
			printf("CUDA Driver Version / Runtime Version    %d.%d / %d.%d\n",
				driverVersion / 1000, (driverVersion % 100) / 10,
				runtimeVersion / 1000, (runtimeVersion % 100) / 10);
			printf("CUDA Capability Major/Minor version number:   %d.%d\n",
				deviceProp.major, deviceProp.minor);
			printf("\n");

			int size = deviceProp.totalGlobalMem / (1024.0 * 1024.0 * 1024.0);
			printf("GPU Memory:   %d G\n",
				size);
			printf("\n");

			gpu_mem = deviceProp.totalGlobalMem / (1024.0 * 1024.0 * 1024.0);
		}
	}

	return true;

}