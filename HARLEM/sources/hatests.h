/*! file hatests.h

    Axxiliary classes to test functionality of HARLEM 

    \author Igor Kurnikov 
    \date 2006

*/
#if !defined(HATESTS_H)
#define HATESTS_H


class HaTests
{
public:
	static void calc_polar_gcontr();
	static void save_grp_oper_mat();
	static void calc_polar_contr_f();
	static void read_polar_contr();
	static void calc_beta_contr_2idx();
	static void read_beta_contr_2idx();
	static void calc_polar_contr_2idx();
	static void read_polar_contr_2idx();
	static void test_oper_1();
	static void test_oper_2();
	static void test_qcmod_1();
	static void dump_mol_info();
	static void dump_gauss_bcommon();
	static void dump_overlap();
    static void dump_overlap2();
	static void test_min_1();
	static void test_graph_1();
	static void test_python_1();
	static void model_mc_calc();
};

#endif

