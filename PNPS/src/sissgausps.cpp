//
// C++ Implementation: poissonsolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "sissgausps.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"

#include <stdlib.h>
#include <limits.h>

#define HALTON_MAX_DIMENSION 1229
static const int prime_numbers[HALTON_MAX_DIMENSION] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
 661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
 8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
 8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
 8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
 8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
 9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
 9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
 9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973
};

SISSGausPS::SISSGausPS()
 : GenericSolver()
{
	InitZero();
}


SISSGausPS::~SISSGausPS()
{
	Clear();
}
int SISSGausPS::InitZero()
{
	World=NULL;
	/*FILE *in=fopen("primenum.dat","r");
	fscanf(in,"%d",&NumOfPrimeNumbers);
	pnpPrint("Load %d prime numbers\n",NumOfPrimeNumbers);
	int i;
	PrimeNumbers=new int[NumOfPrimeNumbers];
	for(i=0;i<NumOfPrimeNumbers;i++)
	{
		fscanf(in,"%d",&(PrimeNumbers[i]));
	}
	fclose(in);*/
	return EXIT_SUCCESS;
}
int SISSGausPS::Clear()
{
	DeleteCVecArray(DielMat,lGS_XYZ);
	return EXIT_SUCCESS;
}
int SISSGausPS::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	GridScale = World->GridScale;
	lGS_X = World->GridSize[0]-2;
	lGS_Y = World->GridSize[1]-2;
	lGS_Z = World->GridSize[2]-2;
	lGS_XY = lGS_X*lGS_Y;
	lGS_XYZ = lGS_XY*lGS_Z;
	DMepslen=2*lGS_XY+1;
	
	
	GS_X = World->GridSize[0];
	GS_Y = World->GridSize[1];
	GS_Z = World->GridSize[2];
	GS_XY = GS_X*GS_Y;
	GS_XYZ = GS_XY*GS_Z;
	
	if(World->Potential == NULL){
		if(!(World->Potential = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
	return EXIT_SUCCESS;
}
int SISSGausPS::ShowParameters()
{
	DbgPrint2("SISSGausPS::ShowParameters\n");
	return EXIT_SUCCESS;
}
int SISSGausPS::ShowProperties()
{
	DbgPrint2("SISSGausPS::ShowProperties\n");
	return EXIT_SUCCESS;
}
int SISSGausPS::InitSolver()
{
	DbgPrint2("SISSGausPS::InitSolver\n");
	float *Q=World->NIndexing->Q;
	float *Eps=World->NIndexing->Eps;
	unsigned int specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex* NIndex=World->NIndexing->NIndex;
	
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	
	int dielectricXS,dielectricYS,dielectricZS,dielectricZSSUM;
	int dielectricXmS,dielectricYmS,dielectricZmS;
	float QstS;
	int i,count;
	int ix,iy,iz;
	int lGrdPnt,GrdPnt;
	int lEps[NodeIndexMaxValues];
	DbgPrint2("Making diel.const. integer:\n");
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		lEps[i]=-int(Eps[i]*EPKT+0.5);
		DbgPrint2("lEps[%d]=%d Eps[%d]=%f\n",i,-lEps[i],i,Eps[i]*EPKT);
	}
	//lEps[0]=1;
	//Get Qst
	Qnum=World->NIndexing->QNum;
	DbgPrint2("Number of static  charges: %d\n",Qnum);
	Qpos=new int[Qnum];
	Qval=new float[Qnum];
	QMult=new DielMatElt*[lGS_XYZ];
	for(lGrdPnt=0;lGrdPnt<lGS_XYZ;lGrdPnt++)
	{
		QMult[lGrdPnt]=new DielMatElt[Qnum];
		for(i=0;i<Qnum;i++)
			QMult[lGrdPnt][i]=0;
	}
	count=0;
	for(ix=1;ix<GS_X-1;ix++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(iz=1;iz<GS_Z-1;iz++)
	{
		GrdPnt = ix+iy*GS_X+iz*GS_XY;
		lGrdPnt=ix-1+(iy-1)*lGS_X+(iz-1)*lGS_XY;
		if((NIndex[GrdPnt]&specChargeMask)==specChargeMask)
		{
			Qpos[count]=lGrdPnt;
			Qval[count]=Q[count]/EPKT;
			QMult[lGrdPnt][count]=1;
		}
	}
	//fill matrix
	DielMat=new DielMatElt*[lGS_XYZ];
	for(lGrdPnt=0;lGrdPnt<lGS_XYZ;lGrdPnt++)
	{
		DielMat[lGrdPnt]=new DielMatElt[DMepslen];
		for(i=0;i<DMepslen;i++)
			DielMat[lGrdPnt][i]=0;
	}
	rDielMat=new DielMatElt*[lGS_XYZ];
	for(lGrdPnt=0;lGrdPnt<lGS_XYZ;lGrdPnt++)
	{
		rDielMat[lGrdPnt]=new DielMatElt[lGS_XYZ];
		for(i=0;i<lGS_XYZ;i++)
			rDielMat[lGrdPnt][i]=0;
	}
	int dmMGSXY=0;
	int dmMGSX=lGS_XY-lGS_X;
	int dmMO=lGS_XY-1;
	int dmCore=lGS_XY;
	int dmPO=lGS_XY+1;
	int dmPGSX=lGS_XY+lGS_X;
	int dmPGSXY=2*lGS_XY;
	//fill internal matrix
	for(ix=1;ix<lGS_X-1;ix++)
		for(iy=1;iy<lGS_Y-1;iy++)
			for(iz=1;iz<lGS_Z-1;iz++)
	{
		lGrdPnt=ix+iy*lGS_X+iz*lGS_XY;
		GrdPnt=ix+1+(iy+1)*GS_X+(iz+1)*GS_XY;
		DielMat[lGrdPnt][dmPO]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		DielMat[lGrdPnt][dmMO]=lEps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		DielMat[lGrdPnt][dmPGSX]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		DielMat[lGrdPnt][dmMGSX]=lEps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		DielMat[lGrdPnt][dmPGSXY]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		DielMat[lGrdPnt][dmMGSXY]=lEps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		DielMat[lGrdPnt][dmCore]=DielMat[lGrdPnt][dmPO]+DielMat[lGrdPnt][dmMO]+DielMat[lGrdPnt][dmPGSX]+DielMat[lGrdPnt][dmMGSX]+DielMat[lGrdPnt][dmPGSXY]+DielMat[lGrdPnt][dmMGSXY];
		DielMat[lGrdPnt][dmCore]=-DielMat[lGrdPnt][dmCore];
		
		rDielMat[lGrdPnt][lGrdPnt+1]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		rDielMat[lGrdPnt][lGrdPnt-1]=lEps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		rDielMat[lGrdPnt][lGrdPnt+lGS_X]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		rDielMat[lGrdPnt][lGrdPnt-lGS_X]=lEps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		rDielMat[lGrdPnt][lGrdPnt+lGS_XY]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		rDielMat[lGrdPnt][lGrdPnt-lGS_XY]=lEps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		rDielMat[lGrdPnt][lGrdPnt]=rDielMat[lGrdPnt][lGrdPnt+1] + rDielMat[lGrdPnt][lGrdPnt-1] + rDielMat[lGrdPnt][lGrdPnt+lGS_X] + rDielMat[lGrdPnt][lGrdPnt-lGS_X] + rDielMat[lGrdPnt][lGrdPnt+lGS_XY] + rDielMat[lGrdPnt][lGrdPnt-lGS_XY];
		rDielMat[lGrdPnt][lGrdPnt]=-rDielMat[lGrdPnt][lGrdPnt];
	}
	//Fill Borders, optimize it later
	int sum;
	for(ix=0;ix<lGS_X;ix++)
		for(iy=0;iy<lGS_Y;iy++)
			for(iz=0;iz<lGS_Z;iz++)
	{
		lGrdPnt=ix+iy*lGS_X+iz*lGS_XY;
		GrdPnt=ix+1+(iy+1)*GS_X+(iz+1)*GS_XY;
		if(ix!=lGS_X-1)
			DielMat[lGrdPnt][dmPO]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		if(ix!=0)
			DielMat[lGrdPnt][dmMO]=lEps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		if(iy!=lGS_Y-1)
			DielMat[lGrdPnt][dmPGSX]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		if(iy!=0)
			DielMat[lGrdPnt][dmMGSX]=lEps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		if(iz!=lGS_Z-1)
			DielMat[lGrdPnt][dmPGSXY]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		if(iz!=0)
			DielMat[lGrdPnt][dmMGSXY]=lEps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		DielMat[lGrdPnt][dmCore]=DielMat[lGrdPnt][dmPO]+DielMat[lGrdPnt][dmMO]+DielMat[lGrdPnt][dmPGSX]+DielMat[lGrdPnt][dmMGSX]+DielMat[lGrdPnt][dmPGSXY]+DielMat[lGrdPnt][dmMGSXY];
		DielMat[lGrdPnt][dmCore]=-DielMat[lGrdPnt][dmCore];
		
		
		rDielMat[lGrdPnt][lGrdPnt]=0;
		if(ix!=lGS_X-1)
		{
			rDielMat[lGrdPnt][lGrdPnt+1]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt+1];
		}
		if(ix!=0)
		{
			rDielMat[lGrdPnt][lGrdPnt-1]=lEps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt-1];
		}
		if(iy!=lGS_Y-1)
		{
			rDielMat[lGrdPnt][lGrdPnt+lGS_X]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt+lGS_X];
		}
		if(iy!=0)
		{
			rDielMat[lGrdPnt][lGrdPnt-lGS_X]=lEps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt-lGS_X];
		}
		if(iz!=lGS_Z-1)
		{
			rDielMat[lGrdPnt][lGrdPnt+lGS_XY]=lEps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt+lGS_XY];
		}
		if(iz!=0)
		{
			rDielMat[lGrdPnt][lGrdPnt-lGS_XY]=lEps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
			rDielMat[lGrdPnt][lGrdPnt]+=rDielMat[lGrdPnt][lGrdPnt-lGS_XY];
		}
		rDielMat[lGrdPnt][lGrdPnt]=-rDielMat[lGrdPnt][lGrdPnt];
	}
	
	/*PrintDielMatCurState();
	for(lGrdPnt=0;lGrdPnt<lGS_XYZ;lGrdPnt++)
	{
		FindAndDevByComMult(lGrdPnt);
	}
	pnpPrint("after FindAndDevByComMult\n");*/
	//PrintDielMatCurState();
	PrintRDielMatCurState();
	
	return EXIT_SUCCESS;
}
int SISSGausPS::PrintDielMat()
{
	int i;
	int ilGrdPnt,jlGrdPnt;
	int dmMGSXY=0;
	int dmMGSX=lGS_XY-lGS_X;
	int dmMO=lGS_XY-1;
	int dmCore=lGS_XY;
	int dmPO=GS_XY+1;
	int dmPGSX=lGS_XY+lGS_X;
	int dmPGSXY=2*lGS_XY;
	int ival;
	pnpPrint("DielMat, as store:\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		for(i=0;i<DMepslen;i++)
		{
			ival=DielMat[ilGrdPnt][i];
			pnpPrint("%4d",ival);
			if(i==DMepslen-1)pnpPrint("\n");
			else pnpPrint(" ");
		}
	}
	pnpPrint("\n");
	pnpPrint("DielMat, Matrix view:\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		int Left=ilGrdPnt-lGS_XY;
		int Right=ilGrdPnt+lGS_XY;
		for(jlGrdPnt=0;jlGrdPnt<lGS_XYZ;jlGrdPnt++)
		{
			ival=0;
			i=jlGrdPnt-ilGrdPnt+lGS_XY;
			if(i>=0&&i<DMepslen)ival=DielMat[ilGrdPnt][i];
			//pnpPrint("%4d,%2d",ival,i);
			pnpPrint("%4d",ival);
			if(jlGrdPnt==lGS_XYZ-1)pnpPrint("\n");
			else pnpPrint(" ");
			
		}
	}
	pnpPrint("\n");
	pnpPrint("Qst:\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		bool hasvalue=false;
		for(i=0;i<Qnum;i++)
		{
			if(QMult[ilGrdPnt][i]!=0)
			{
				if(hasvalue)
					pnpPrint("+");
				else
					hasvalue=true;
				pnpPrint("%d*(Q_%d=%f)",QMult[ilGrdPnt][i],Qpos[i],Qval[i]);
			}
		}
		if(!hasvalue)
			pnpPrint("0\n");
	}
	pnpPrint("\n");
	return EXIT_SUCCESS;
}

int SISSGausPS::PrintDielMatCurState()
{
	int i;
	int ilGrdPnt,jlGrdPnt;
	int ival;
	pnpPrint("|DielMat||Qst|\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		int Left=ilGrdPnt-lGS_XY;
		int Right=ilGrdPnt+lGS_XY;
		pnpPrint("%4d | ",ilGrdPnt);
		
		for(jlGrdPnt=0;jlGrdPnt<lGS_XYZ;jlGrdPnt++)
		{
			ival=0;
			i=jlGrdPnt-ilGrdPnt+lGS_XY;
			if(i>=0&&i<DMepslen)ival=DielMat[ilGrdPnt][i];
			//pnpPrint("%4d,%2d",ival,i);
			pnpPrint("%4d ",ival);
		}
		pnpPrint("| ",ival);
		bool hasvalue=false;
		for(i=0;i<Qnum;i++)
		{
			if(QMult[ilGrdPnt][i]!=0)
			{
				if(hasvalue)
					pnpPrint("+");
				else
					hasvalue=true;
				pnpPrint("%d*(Q_%d=%f)",QMult[ilGrdPnt][i],Qpos[i],Qval[i]);
			}
		}
		if(!hasvalue)
			pnpPrint("0");
		pnpPrint("\n");
	}
	pnpPrint("\n");
	return EXIT_SUCCESS;
}
int SISSGausPS::PrintRDielMat()
{
	int i;
	int ilGrdPnt,jlGrdPnt;
	int ival;
	pnpPrint("|rDielMat|\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		for(jlGrdPnt=0;jlGrdPnt<lGS_XYZ;jlGrdPnt++)
		{
			pnpPrint("%d",rDielMat[ilGrdPnt][jlGrdPnt]);
			if(jlGrdPnt!=lGS_XYZ-1)
				pnpPrint("\t");
			else
				pnpPrint("\n");
		}
	}
	pnpPrint("\n");
	return EXIT_SUCCESS;
}
int SISSGausPS::PrintRDielMatCurState()
{
	int i;
	int ilGrdPnt,jlGrdPnt;
	int ival;
	pnpPrint("|rDielMat||Qst|\n");
	for(ilGrdPnt=0;ilGrdPnt<lGS_XYZ;ilGrdPnt++)
	{
		pnpPrint("%4d | ",ilGrdPnt);
		
		for(jlGrdPnt=0;jlGrdPnt<lGS_XYZ;jlGrdPnt++)
		{
			pnpPrint("%4d ",rDielMat[ilGrdPnt][jlGrdPnt]);
		}
		pnpPrint("| ",ival);
		bool hasvalue=false;
		for(i=0;i<Qnum;i++)
		{
			if(QMult[ilGrdPnt][i]!=0)
			{
				if(hasvalue)
					pnpPrint("+");
				else
					hasvalue=true;
				pnpPrint("%d*(Q_%d=%f)",QMult[ilGrdPnt][i],Qpos[i],Qval[i]);
			}
		}
		if(!hasvalue)
			pnpPrint("0");
		pnpPrint("\n");
	}
	pnpPrint("\n");
	return EXIT_SUCCESS;
}
int SISSGausPS::AddRowItoRowJGoDown(int I,DielMatElt multI,int J,DielMatElt multJ)
{
	// always dIJ>0and dIJ>=DMepslen
	int i;
	int dIJ=J-I;
	//DielMat
	//MultI
	for(i=0;i<DMepslen;i++)
		DielMat[J][i]*=multJ;
	
	for(i=0;i<DMepslen-dIJ;i++)
	{
		DielMat[J][i]+=multI*DielMat[I][i+dIJ];
	}
	//QMult
	for(i=0;i<Qnum;i++)
	{
		QMult[J][i]=multJ*QMult[J][i]+multI*QMult[I][i];
	}
	return EXIT_SUCCESS;
}
int SISSGausPS::AddRowItoRowJGoUp(int I,DielMatElt multI,int J,DielMatElt multJ)
{
	// always dIJ>0and dIJ>=DMepslen
	int i;
	int dIJ=J-I;
	//DielMat
	//MultI
	for(i=0;i<DMepslen;i++)
		DielMat[J][i]*=multJ;
	
	for(i=-dIJ;i<DMepslen;i++)
	{
		DielMat[J][i]+=multI*DielMat[I][i+dIJ];
	}
	
	//for(i=0;i<DMepslen-dIJ;i++)
	//{
	//i=lGS_XY-dIJ;
	//	DielMat[J][i]+=multI*DielMat[I][i+dIJ];
	//}
	//QMult
	for(i=0;i<Qnum;i++)
	{
		QMult[J][i]=multJ*QMult[J][i]+multI*QMult[I][i];
	}
	return EXIT_SUCCESS;
}
bool SISSGausPS::IsItComMultForRowJ(int J,DielMatElt Mult)
{
	int i;
	for(i=0;i<Qnum;i++)
	{
		if(QMult[J][i]%Mult!=0)
			return false;
	}
	for(i=0;i<DMepslen;i++)
	{
		if(DielMat[J][i]%Mult!=0)
			return false;
	}
	return true;
}
bool SISSGausPS::RIsItComMultForRowJ(int J,DielMatElt Mult)
{
	int i;
	for(i=0;i<Qnum;i++)
	{
		if(QMult[J][i]%Mult!=0)
			return false;
	}
	for(i=0;i<lGS_XYZ;i++)
	{
		if(DielMat[J][i]%Mult!=0)
			return false;
	}
	return true;
}
int SISSGausPS::DevRowJByComMult(int J,DielMatElt Mult)
{
	int i;
	for(i=0;i<DMepslen;i++)
	{
		DielMat[J][i]/=Mult;
	}
	for(i=0;i<Qnum;i++)
	{
		QMult[J][i]/=Mult;
	}
	return EXIT_SUCCESS;
}
int SISSGausPS::RDevRowJByComMult(int J,DielMatElt Mult)
{
	int i;
	for(i=0;i<lGS_XYZ;i++)
	{
		rDielMat[J][i]/=Mult;
	}
	for(i=0;i<Qnum;i++)
	{
		QMult[J][i]/=Mult;
	}
	return EXIT_SUCCESS;
}
int SISSGausPS::FindAndDevByComMult(int J)
{
	int i,val;
	int iMin=-1,iMax=-1;
	for(i=0;i<DMepslen;i++)
	{
		val=abs(DielMat[J][i]);
		
		if(val>0&&iMin==-1)
		{
			iMin=val;
			iMax=val;
		}
		if(val<iMin)iMin=val;
		if(val>iMax)iMax=val;
	}
	int num=0;
	while(num<NumOfPrimeNumbers && PrimeNumbers[num]<=iMax)
	{
		while(IsItComMultForRowJ(J,PrimeNumbers[num]))
			DevRowJByComMult(J,PrimeNumbers[num]);
		num++;
	}
	 
	return EXIT_SUCCESS;
}
int SISSGausPS::RFindAndDevByComMult(int J)
{
	int i,val;
	int iMin=-1,iMax=-1;
	for(i=0;i<lGS_XYZ;i++)
	{
		val=abs(DielMat[J][i]);
		
		if(val>0&&iMin==-1)
		{
			iMin=val;
			iMax=val;
		}
		if(val<iMin)iMin=val;
		if(val>iMax)iMax=val;
	}
	int num=0;
	/*while(num<NumOfPrimeNumbers and PrimeNumbers[num]<=iMax)
	{
		while(IsItComMultForRowJ(J,PrimeNumbers[num]))
			DevRowJByComMult(J,PrimeNumbers[num]);
		num++;
	}*/
	while(num<HALTON_MAX_DIMENSION && prime_numbers[num]<=iMax)
	{
		while(RIsItComMultForRowJ(J,prime_numbers[num]))
			RDevRowJByComMult(J,prime_numbers[num]);
		num++;
	}
	 
	return EXIT_SUCCESS;
}
int SISSGausPS::RexchangeRows(int I, int J)
{
	int i;
	int ilGrdPnt,jlGrdPnt;
	int ival;
	
	DielMatElt *tmp;
	
	tmp=rDielMat[I];
	rDielMat[I]=rDielMat[J];
	rDielMat[J]=tmp;
	
	tmp=QMult[I];
	QMult[I]=QMult[J];
	QMult[J]=tmp;
	
	return EXIT_SUCCESS;
}
int SISSGausPS::RAddRowItoRowJ(int I,DielMatElt multI,int J,DielMatElt multJ)
{
	// always dIJ>0and dIJ>=DMepslen
	int i;
	//DielMat
	for(i=0;i<lGS_XYZ;i++)
	{
		rDielMat[J][i]=multJ*rDielMat[J][i]+multI*rDielMat[I][i];
	}
	//QMult
	for(i=0;i<Qnum;i++)
	{
		QMult[J][i]=multJ*QMult[J][i]+multI*QMult[I][i];
	}
	return EXIT_SUCCESS;
}
int SISSGausPS::Solve()
{
	int i;
	int I,J;
	int K1,K2;
	
	DbgPrint2("SISSGausPS::Solve\n");
	PrintRDielMat();
	//prepare matrix
	DielMatElt **rowstmp=new DielMatElt*[lGS_XY];
	for(I=0;I<lGS_XY;I++)
	{
		rowstmp[I]=rDielMat[I];
	}
	for(I=0;I<lGS_XYZ-lGS_XY;I++)
	{
		rDielMat[I]=rDielMat[I+lGS_XY];
	}
	for(I=0;I<lGS_XY;I++)
	{
		rDielMat[I+lGS_XYZ-lGS_XY]=rowstmp[I];
	}
	PrintRDielMat();
	
	for(I=lGS_XYZ-lGS_XY;I<lGS_XYZ;I++)
	{
		pnpPrint(" I=%d\n",I);
		for(J=0;J<lGS_XYZ-lGS_XY;J++)
		{
			if(rDielMat[I][J]!=0)
			{
				pnpPrint(" J=%d; %d;",J,rDielMat[J][J]);
				if(rDielMat[I][J]%rDielMat[J][J]==0)
				{
					pnpPrint(" J is multiply of I\n");
					RAddRowItoRowJ(J,-rDielMat[I][J]/rDielMat[J][J],I,1);
				}
				else if(rDielMat[J][J]%rDielMat[I][J]==0)
				{
					//pnpPrint(" I is multiply of J\n");
					RAddRowItoRowJ(J,1,I,-rDielMat[J][J]/rDielMat[I][J]);
				}
				else
				{
					//pnpPrint("heavy case\n");
					RAddRowItoRowJ(J,-rDielMat[I][J],I,rDielMat[J][J]);
				}
			}
		}
		RFindAndDevByComMult(I);
		//PrintRDielMat();
	}
	PrintRDielMat();
	return EXIT_SUCCESS;
	//optimized last GSXY
	//\\\
	pnpPrint("optimized last GSXY\n");
	for(I=lGS_XYZ-lGS_XY;I<lGS_XYZ;I++)
	{
		pnpPrint(" I=%d\n",I);
		//order
		for(K1=I;K1<lGS_XYZ;K1++)
		{
			for(K2=K1+1;K2<lGS_XYZ;K2++)
			{
				if(rDielMat[K1][I]==0&&rDielMat[K2][I]!=0)
				{
					RexchangeRows(K1,K2);
				}
				else if(abs(rDielMat[K2][I])<abs(rDielMat[K1][I])&&rDielMat[K2][I]>0)
				{
					RexchangeRows(K1,K2);
				}
			}
		}
		PrintRDielMat();
		for(K1=lGS_XYZ-1;K1>I;K1--)
		{
			if(rDielMat[K1][I]!=0)
			{
				for(K2=K1-1;K2>=I;K2--)
					if(rDielMat[K1][I]%rDielMat[K2][I]==0)
						RAddRowItoRowJ(K2,-rDielMat[K1][I]/rDielMat[K2][I],K1,1);
			}
			if(rDielMat[K1][I]!=0)
			{
				RAddRowItoRowJ(I,-rDielMat[K1][I],K1,rDielMat[I][I]);
			}
		}
		RFindAndDevByComMult(I);
		PrintRDielMat();
	}
	PrintRDielMat();
	//PrintDielMatCurState();
	return EXIT_SUCCESS;
	
	DbgPrint2("SISSGausPS::Solve\n");
	PrintDielMatCurState();
	//go down
	
	
	int dmMGSXY=0;
	int dmMGSX=lGS_XY-lGS_X;
	int dmMO=lGS_XY-1;
	int dmCore=lGS_XY;
	int dmPO=lGS_XY+1;
	int dmPGSX=lGS_XY+lGS_X;
	int dmPGSXY=2*lGS_XY;
	//lGS_XYZ
	//pnpPrint("int min max %d %d\n",INT_MIN,INT_MAX);return EXIT_FAILURE;
	//-2147483648 2147483647
	//go down
	for(I=0;I<lGS_XYZ;I++)
	{
		if(DielMat[I][dmCore]==0)
		{
			pnpPrint("Row %d: Ups core element is 0, can not do it yet\n",I);
			PrintDielMatCurState();
			return EXIT_FAILURE;
		}
		
		for(J=I+1;J<I+lGS_XY+1&&J<lGS_XYZ;J++)
		{
			int dIJ=J-I;
			if(DielMat[J][dmCore-dIJ]!=0)
			{
				//pnpPrint("%d",J);
				//J is multiply of I
				if(DielMat[J][dmCore-dIJ]%DielMat[I][dmCore]==0)
				{
					//pnpPrint(" J is multiply of I\n");
					AddRowItoRowJGoDown(I,-DielMat[J][dmCore-dIJ]/DielMat[I][dmCore],J,1);
				}
				else if(DielMat[I][dmCore]%DielMat[J][dmCore-dIJ]==0)
				{
					//pnpPrint(" I is multiply of J\n");
					AddRowItoRowJGoDown(I,1,J,-DielMat[I][dmCore]/DielMat[J][dmCore-dIJ]);
				}
				else
				{
					//pnpPrint("heavy case\n");
					AddRowItoRowJGoDown(I,-DielMat[J][dmCore-dIJ],J,DielMat[I][dmCore]);
				}
				FindAndDevByComMult(J);
				//pnpPrint("\n");
			}
		}
	}
	PrintDielMatCurState();
	//Go Up
	for(I=lGS_XYZ-1;I>=0;I--)
	{
		if(DielMat[I][dmCore]==0)
		{
			pnpPrint("Row %d: Ups core element is 0, can not do it yet\n",I);
			PrintDielMatCurState();
			return EXIT_FAILURE;
		}
		pnpPrint("I=%d\n",I);
		for(J=I-1;J>I-lGS_XY-1&&J>=0;J--)
		{
			int dIJ=J-I;
			if(DielMat[J][dmCore-dIJ]!=0)
			{
				pnpPrint("%d",J);
				//J is multiply of I
				if(DielMat[J][dmCore-dIJ]%DielMat[I][dmCore]==0)
				{
					pnpPrint(" J is multiply of I\n");
					AddRowItoRowJGoUp(I,-DielMat[J][dmCore-dIJ]/DielMat[I][dmCore],J,1);
				}
				else if(DielMat[I][dmCore]%DielMat[J][dmCore-dIJ]==0)
				{
					pnpPrint(" I is multiply of J\n");
					AddRowItoRowJGoUp(I,1,J,-DielMat[I][dmCore]/DielMat[J][dmCore-dIJ]);
				}
				else
				{
					pnpPrint("heavy case\n");
					AddRowItoRowJGoUp(I,-DielMat[J][dmCore-dIJ],J,DielMat[I][dmCore]);
				}
				FindAndDevByComMult(J);
				pnpPrint("\n");
			}
		}
		PrintDielMatCurState();
	}
	PrintDielMatCurState();
	pnpPrint("int min max %d %d\n",INT_MIN,INT_MAX);
	//Go Up
	
	//AddRowItoRowJGoDown(0,1,1,1);
	
	//AddRowItoRowJGoDown(0,1,3,1);
	//AddRowItoRowJGoDown(0,1,9,1);
	//PrintDielMatCurState();
	return EXIT_FAILURE;
}

