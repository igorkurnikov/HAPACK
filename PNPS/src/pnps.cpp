/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2005                                      *
 ***************************************************************************/

#include "pnpdebug.h"

#include <omp.h>

#include <string>
using namespace std;

#include "Python.h"

extern "C" {
	extern void init_llpnps(void);
}


#if 0
// !defined(HARLEM_MOD)
int main(int argc, char *argv[])
{
	int argc4py = 0;
	char **argv4py = NULL;
	int threads_number=-1;
	bool printUsage = false;
	bool interactive_mode = false;
	Py_Initialize();
	init_llpnps();

	if (argc > 1)
	{
		int argc4py = 0;
		char **argv4py = new char*[argc];

		int num_threads = 1;

		for (int i = 0; i<argc; i++) {
			string arg(argv[i]);

			//strip off arguments which should not make to python
			if (arg == "-nt")
			{

				if ((i + 1 >= argc) || (!isdigit(argv[i + 1][0])))
				{
					argc4py = 0;
					printUsage = true;//print usage and exit
					printf("Error: Incorrect command line format\n");
					break;
				}
				++i;

				num_threads = atoi(argv[i]);

				//set openmp threads
#ifdef _OPENMP
				printf("Set threads number to %d\n\n", num_threads);
				omp_set_num_threads(num_threads);
#else
				if (num_threads>1)
					printf("Warning: This version was compiled without openmp support!\n");

				num_threads = 1;
				printf("Set threads number to %d\n\n", num_threads);
#endif
				continue;
			}
			if (arg == "-i")
			{
				interactive_mode = true;
				continue;
			}
			argv4py[argc4py] = argv[i];
			++argc4py;
	    }

#ifdef _OPENMP
		threads_number = omp_get_max_threads();
		printf("Max threads number is %d\n\n", threads_number);
#else
		threads_number = 1;
#endif

		if (argc4py > 1)
		{
			//execute user script
			Py_Main(argc4py, argv4py);
		}
		else
		{
			//print usage and exit
			//printUsage = true;
			Py_Main(argc4py, argv4py);
		}
		delete[] argv4py;
    }
	else
	{
		//print usage and exit
		printUsage = true;
	}

	if (interactive_mode)
	{
		Py_Main(argc4py, argv4py);
	}
	else if (printUsage)
	{
		printf("Usage:\n");
		printf("\t%s [-nt <number of threads>] [-i] <run_script>]\n", argv[0]);
		printf("\t-i interactive mode\n");
	}

	Py_Finalize();
	return EXIT_SUCCESS;
}
#endif
