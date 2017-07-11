#include <time.h>
#include <sys/time.h>

double get_wall_time()
{
#if 1
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
#endif
}

double get_cpu_time()
{
	return (double)clock() / CLOCKS_PER_SEC;
}

class Stopwatch {
public: 
	Stopwatch() {
		sumWallTime = 0;
		startWallTime = -1;
		lastWallTime = -1;
		
		sumCpuTime = 0;
		startCpuTime = -1;
		lastCpuTime = -1;
		running = false;
	}
	
	void start() {
		startWallTime = get_wall_time();
		startCpuTime = get_cpu_time();
		running = true;
	}
	
	void stop() {
		double endWallTime = get_wall_time();
		sumWallTime = sumWallTime + (endWallTime - startWallTime);
		lastWallTime = endWallTime - startWallTime;
		
		double endCpuTime = get_cpu_time();
		sumCpuTime = sumCpuTime + (endCpuTime - startCpuTime);
		lastCpuTime = endCpuTime - startCpuTime;
		
		running = false;
	}
	
	void getTime(double* times) {
		times[0] = lastWallTime;
		times[1] = lastCpuTime;
	}
	
	void addTime(double* times) {
		times[0] += lastWallTime;
		times[1] += lastCpuTime;
	}
	
	void getTotal(double* times) {
		times[0] = sumWallTime;
		times[1] = sumCpuTime;
	}
	
private:
	double sumWallTime;
	double startWallTime;
	double lastWallTime;
	
	double sumCpuTime;
	double startCpuTime;
	double lastCpuTime;
	bool running;
};
