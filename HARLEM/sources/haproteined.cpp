/*! \file haproteined.cpp

    Classes to perform Essential Dynamics and clustering analysis

    \author Igor Kurnikov  
    \date 2011-

*/  

#define HAPROTEINED_CPP

#include <mpi.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "hastl.h"
#include "haio.h"

#include "harlemapp.h"

#include "hacompmod.h"
#include "haatgroup.h"

#include "hamolmech.h"
#include "mm_traj_anal.h"
#include "haproteined.h"

#include "tinyxml.h"

#include "hamolset.h"

IntStrMap SetCollectCrdAnalTypeLbls()
{
	IntStrMap lmap;
	
	lmap[CollectCrdAnalType::CCRD_SPATIAL_PLATO ] =  "Spatial PLATO";
	lmap[CollectCrdAnalType::CCRD_SPATIAL_RV    ] =  "Spatial RV";
	lmap[CollectCrdAnalType::CCRD_TEMPORAL_PCA  ] =  "Temporal PCA";
	lmap[CollectCrdAnalType::CCRD_TEMPORAL_ICA  ] =  "Temporal ICA";
	lmap[CollectCrdAnalType::CCRD_TEMPORAL_SFA  ] =  "Temporal SFA";
	lmap[CollectCrdAnalType::CCRD_TEMPORAL_ISFA ] =  "Temporal ISFA";
	
	return lmap;
}

IntStrMap CollectCrdAnalType::labels = SetCollectCrdAnalTypeLbls();  

CollectCrdAnalType::CollectCrdAnalType()
{
	v_= CCRD_TEMPORAL_PCA;	
}
	
CollectCrdAnalType::~CollectCrdAnalType()
{

}

int CollectCrdAnalType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int CollectCrdAnalType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}


CollectCrdAnalMod::CollectCrdAnalMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_CLUSTER_ANAL, new_phost_mset)
{
	p_md_traj = new MDTrajectory(new_phost_mset);
	md_traj_fname_plato = "traj.plato";
	plato_input_fname   = "plato_inp.xml";
	plato_output_fname  = "plato_out.xml";
	plato_run_fname     = "plato_run.py";
	plato_root_dir      = "C:/MY_PROJECTS/GEORG_ESSENTIAL_DYNAMICS";
	time_proj_fname     = "time_proj.dat";

	num_clusters  = 4;        
	num_eig_vec   = 4;        
	num_time_proj = 4;        
	sim_matrix_flag = TRUE;

	load_vector_type = LOAD_WVECTORS;
}
	
CollectCrdAnalMod::~CollectCrdAnalMod()
{
	delete p_md_traj;
}

void CollectCrdAnalMod::SetActiveAtomGroup(const std::string& atgrp_name)
{
	active_atgrp = atgrp_name;
}

std::string CollectCrdAnalMod::GetActiveAtomGroupName()
{
	return active_atgrp;
}

AtomGroup* CollectCrdAnalMod::GetActiveAtomGroup()
{
	MolSet* pmset = this->GetMolSet();
	return pmset->GetAtomGroupByID(active_atgrp.c_str());
}

int CollectCrdAnalMod::ConvertMDTrajToPlato()
{
	char buf[256];
	int is_md_opened = FALSE;
	try
	{
		AtomGroup* p_at_grp = GetActiveAtomGroup();
		if( p_at_grp == NULL || p_at_grp->GetNAtoms() == 0 ) throw std::runtime_error( " Atom Group " + active_atgrp + " is non-existent or empty" );

		int ires = p_md_traj->Open();
		if( !ires) throw std::runtime_error(" Error to open MD trajectory: " + p_md_traj->CrdFileName );
		is_md_opened = TRUE;

		std::ofstream plato_s( md_traj_fname_plato.c_str() );
		if( !plato_s ) throw std::runtime_error(" Error to open PLATO trajectory file for writing: " + md_traj_fname_plato );

		int nat = p_at_grp->GetNAtoms();
		int i;
		for( i = 0; i < 3*nat; i++)
		{
			if( i > 0 ) plato_s << ";";
			plato_s << "3";
		}
		plato_s << std::endl;
		
		AtomIteratorAtomGroup aitr(p_at_grp);
		ires = p_md_traj->ReadNextFrame();
		while(ires)
		{
			HaAtom* aptr;
			int first_atom = TRUE;
			for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
			{
				if(!first_atom) plato_s << ";";
				sprintf(buf,"%8.4f",aptr->GetX_Ang());
				plato_s << buf << ";";
				sprintf(buf,"%8.4f",aptr->GetY_Ang());
				plato_s << buf << ";";
				sprintf(buf,"%8.4f",aptr->GetZ_Ang());
				plato_s << buf;
			
				first_atom = FALSE;
			}
			plato_s << std::endl;
			ires = p_md_traj->ReadNextFrame();
		}
	}
	catch ( std::exception& ex )
	{
		if( is_md_opened ) p_md_traj->Close();
		PrintLog(" Error in CollectCrdAnalMod::ConvertMDTrajToPlato() \n");
		PrintLog(" %s \n",ex.what());
		return FALSE;
	}
	PrintLog(" Finish to convert trajectory to PLATO format \n");
	return TRUE;	
}

int CollectCrdAnalMod::LoadPlatoOutputFile()
{
	char buf[256];
	std::string str;
	try
	{
		AtomGroup* p_at_grp = GetActiveAtomGroup();
		if( p_at_grp == NULL ) throw std::runtime_error( "Active Atom Group is not set");
		int nat = p_at_grp->GetNAtoms();

		FILE* fplato_out = fopen( plato_output_fname.c_str(),"r");
		if( !fplato_out ) throw std::runtime_error( "Error to open file " + plato_output_fname);
			
		TiXmlDocument doc;
		bool bres = doc.LoadFile(fplato_out); 
		if(!bres) throw std::runtime_error( "XML file " + plato_output_fname + " is invalid ");

		int ires;
		MolSet* pmset = GetMolSet();
		const TiXmlElement* root_element;
		const TiXmlElement* output_element;

		root_element = doc.FirstChildElement();
		if( root_element == NULL ) throw std::runtime_error( " No root element in XML file " + plato_output_fname );
		
		output_element = root_element->FirstChildElement("output");
		if( output_element == NULL)  throw std::runtime_error( " No element with name output in XML file " + plato_output_fname );
	
		const TiXmlElement* spatial_element = output_element->FirstChildElement("spatial");
		if( spatial_element )
		{
			const TiXmlElement* clusters_element = spatial_element->FirstChildElement("clusters");
			if( clusters_element )
			{
				PrintLog("Load Atom Clusters information \n");
				int nc;
				ires = clusters_element->GetIntAttribute("nclusters",&nc);
				if(ires == EXIT_FAILURE ) throw std::runtime_error( " Error to read the number of clusters " + plato_output_fname );
				std::vector<std::string> str_arr;
				std::string text = clusters_element->GetText();
				boost::split(str_arr,text,boost::is_any_of(","),boost::token_compress_on);
				int nat_pl = str_arr.size();	
				if( nat_pl != nat ) 
				{
					sprintf(buf,"N atoms in clusters element %d is not equal to N atoms in active group ",nat_pl,nat);
					throw std::runtime_error( buf );
				}
				std::vector<AtomGroup*> cluster_groups;
				int ic;
				for(ic = 1; ic <= nc; ic++ )
				{
					sprintf(buf,"CLUSTER_%d_%d",ic,nc);
					std::string cl_grp_name = buf;
					AtomGroup* pgrp = pmset->GetAtomGroupByID(cl_grp_name.c_str());
					if( pgrp == NULL) 
					{
						pgrp = pmset->AddAtomGroup(cl_grp_name.c_str());
					}
					pgrp->clear();
					cluster_groups.push_back(pgrp);
				}

				int i;
				for(i = 0; i < nat; i++)
				{
					ic = atoi(str_arr[i].c_str());
					HaAtom* aptr = (*p_at_grp)[i];
					if( ic == 0 ) if( nat_pl != nat ) throw std::runtime_error(" Error to convert cluster number in cluster element "); 
					if( ic < 0 || ic > nc ) throw std::runtime_error(" Invalid cluster index " + str_arr[i] ); 
					cluster_groups[ic-1]->push_back(aptr);
				}
			}
		}

		const TiXmlElement* temporal_element = output_element->FirstChildElement("temporal");
		if( temporal_element )
		{
			std::vector<std::string> str_arr;
			std::string text;
			const TiXmlElement* eigen_vals_element = temporal_element->FirstChildElement("eigenvalues");
			int neig = 0;
			if( eigen_vals_element )
			{
				PrintLog("Load Eigenvalues \n");
				int nc;
				text = eigen_vals_element->GetText();
				boost::split(str_arr,text,boost::is_any_of(","),boost::token_compress_on);
				neig = str_arr.size();	
				if( neig == 0 )  throw std::runtime_error(" Zero number of eigen values - no temporal info is read from PLATO file "); 
				eigen_vals.clear();
				eigen_vecs.clear();
				int i;
				for(i = 0; i < neig; i++)
				{
					PrintLog(" str_arr[%d] = %s  %12.6f \n",i,str_arr[i].c_str(),atof(str_arr[i].c_str()));
					eigen_vals.push_back( atof(str_arr[i].c_str()) );
				}			
			}

			const TiXmlElement* eigen_vecs_element = temporal_element->FirstChildElement("weightvectors");
			if( eigen_vecs_element )
			{
				PrintLog("Load Eigenvectors \n");
				const TiXmlElement* eigen_vec_element = eigen_vecs_element->FirstChildElement("weightvector");
				while( eigen_vec_element )
				{
					text = eigen_vec_element->GetText();
					boost::split(str_arr,text,boost::is_any_of(","),boost::token_compress_on);
					int nsize_vec = str_arr.size();	
					if( nsize_vec != 3*nat )  
					{
						sprintf(buf," The size of eigenvector array %d is not equal to 3 X number of atoms in the active group %d ",nsize_vec,3*nat);
						throw std::runtime_error(buf);
					}
					HaVec_double evec(3*nat);
					int i;
					double dmax = 0.0;
					for( i = 0; i < nat; i++ )
					{
						double x = atof( str_arr[3*i  ].c_str());
						double y = atof( str_arr[3*i+1].c_str());
						double z = atof( str_arr[3*i+2].c_str());
						evec[3*i  ] = x;
						evec[3*i+1] = y;
						evec[3*i+2] = z;
						double dlen = sqrt( x*x + y*y + z*z);
						if( dlen > dmax ) dmax = dlen;
					}
					if( dmax > 1e-6 )
					{
						for( i = 0; i < nat; i++ )
						{
							evec[3*i  ] /= dmax;
							evec[3*i+1] /= dmax;
							evec[3*i+2] /= dmax;
						}
					}

					eigen_vecs.push_back(evec);
					eigen_vec_element = eigen_vec_element->NextSiblingElement("weightvector");
				}
			} // end weightvectors

			const TiXmlElement* time_projections_element = temporal_element->FirstChildElement("time_projections");
			if( time_projections_element )
			{
				PrintLog("Load Time Projections \n");
				time_projections.clear();
				const TiXmlElement* time_proj_element = time_projections_element->FirstChildElement("time_projection");
				while( time_proj_element )
				{
					text = time_proj_element->GetText();
					boost::split(str_arr,text,boost::is_any_of(","),boost::token_compress_on);
					int nsize = str_arr.size();	
					
					HaVec_double proj_vec(nsize);
					int i;
					for( i = 0; i < nsize; i++ )
					{
						proj_vec[i] = atof( str_arr[i].c_str());
					}
					time_projections.push_back(proj_vec);
					time_proj_element = time_proj_element->NextSiblingElement("time_projection");
				}
				int nproj = time_projections.size();
				PrintLog("The number of time projections read is %d \n", nproj);
				
				if( nproj > 0 ) SaveTimeProjFile();
			} // end time_projections
		}
	}
	catch(std::exception& ex) 
	{
		PrintLog("Error in CollectCrdAnalMod::LoadPlatoOutputFile() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}


int CollectCrdAnalMod::SavePlatoInputFiles()
{
	PrintLog(" Save Input files for PLATO statistical program \n");
	try
	{
		AtomGroup* p_at_grp = GetActiveAtomGroup();
		if( p_at_grp == NULL ) throw std::runtime_error(" No Active Atom Group is set " );

		if( p_at_grp->GetNAtoms() == 0 ) throw std::runtime_error(" Active Atom Group is empty() ");

		std::ofstream os(plato_input_fname.c_str());
		if( os.fail() ) throw std::runtime_error(" Error to open PLATO input file " +  plato_input_fname );
 
		os << "<?xml version=\"1.0\" ?>" << std::endl;
		os << "<trajectory name=\"PLATO trajectory\">" << std::endl;
		os << "  <file dec=\".\" ref=\"" << md_traj_fname_plato << "\" sep=\";\" type=\"text/plato\"/>" << std::endl;
		
		os << "  <description>"  << std::endl;
		
		os << "    <comments>" << std::endl;
		os << "      <comment author=\"HARLEM\" date=\"Today\">" << std::endl;
		os << "        HARLEM generated PLATO input file " << std::endl;
		os << "      </comment>" << std::endl;
		os << "    </comments>" << endl;

		int i,na;
		na = p_at_grp->GetNAtoms();

		os << "    <properties Ndim=\"3\" natoms=\"" << na << "\" nobs=\"99999\" time_step=\"1\" time_unit=\"s\"/>" << std::endl;
		os << "    <labels>";
		for( i = 0; i < na; i++)
		{
			if( ((i % 10) == 0) ) os << std::endl << "      ";
			os << (*p_at_grp)[i]->GetRef();
			if( i == (na - 1) ) os << std::endl;
			else os << ", ";
		}
		os << "    </labels>" << std::endl; 
		os << "  </description>" << std::endl;

		os << "  <input program=\"HARLEM\">" << std::endl;
		os << "    <data>" << std::endl;
		os << "      <transformation>"  << std::endl;
		os << "        none" << std::endl; 
		os << "      </transformation>" << std::endl;
		os << "      <window delta_t=\"1\" end=\"T\" start=\"1\"/>" << std::endl;
		os << "      <subselection nselected=\"" << na << "\">";
		for( i = 0; i < na; i++)
		{
			if( ((i % 10) == 0) ) os << std::endl << "      ";
			os << (i+1);
			if( i == (na - 1) ) os << std::endl;
			else os << ", ";
		}
		os << "      </subselection>" << std::endl;
		os << "    </data>" << std::endl;

		std::string anal_type_str;
		std::string method_str;
		if( anal_type == anal_type.CCRD_SPATIAL_PLATO || anal_type == anal_type.CCRD_SPATIAL_RV )
		{
			anal_type_str = "spatial";
			method_str    = "PLATO";
		}
		else if( anal_type == anal_type.CCRD_SPATIAL_RV )
		{
			anal_type_str = "spatial";
			method_str    = "RV";
		}
		else if( anal_type == anal_type.CCRD_TEMPORAL_ICA )
		{
			anal_type_str = "temporal";
			method_str    = "ICA";
		}
		else if( anal_type == anal_type.CCRD_TEMPORAL_PCA )
		{
			anal_type_str = "temporal";
			method_str    = "PCA";
		}
		else if( anal_type == anal_type.CCRD_TEMPORAL_ISFA )
		{
			anal_type_str = "temporal";
			method_str    = "ISFA";
		}
		else if( anal_type == anal_type.CCRD_TEMPORAL_SFA)
		{
			anal_type_str = "temporal";
			method_str    = "SFA";
		}

		os << "    <analysis type= \"" << anal_type_str << "\">" << std::endl;
		if( anal_type_str == "temporal" )
		{
			os << "      <method TimeFrames=\"1\" nclusters=\"" << num_clusters << "\" ncomp=\"" << num_eig_vec << "\" nneighbors=\"500\" ";
		    os << " nonlinearity=\"none\" nonlinearity_degree=\"1\" preprojection_dim=\"" << num_time_proj << "\">" << std::endl;
		}
		else if( anal_type_str == "spatial" )
		{
			os << "      <method nclusters=\"" << num_clusters << "\"" << " neigenvectors=\"" << num_eig_vec << "\" ncomp=\"" << num_clusters << "\" >" << std::endl;
		}
		os << "        " << method_str << std::endl;
		os << "      </method>" << std::endl;
		os << "      <output_params>" << std::endl;
		os << "        <output_file ref=\"" << plato_output_fname << "\" type=\"text/xml\"/>" << std::endl;
		if( anal_type_str == "temporal" )
		{
		    os << "        <return eigenvalues=\"" << num_eig_vec << "\" time_projections=\"" << num_time_proj << "\" "; 
		    os << " weightvectors=\"" << num_eig_vec << "\"/>" << std::endl;
		}
		else if( anal_type_str == "spatial" )
		{
		    os << "        <return clusters=\"True\" eigenvalues=\"0\" eigenvectors=\"0\" ";
		    os << "similarity_matrix=\"" << (sim_matrix_flag ? "True" : "False") << "\" />" << std::endl;
		}
		os << "      </output_params>" << std::endl;
		os << "    </analysis>" << std::endl;
		os << "  </input>" << std::endl;

		os << "</trajectory>" << std::endl;

		std::ofstream os_run(plato_run_fname.c_str());
		if( os_run.fail() ) throw std::runtime_error(" Error to open PLATO run file " +  plato_run_fname );

		os_run << "import os " << std::endl << std::endl;
		os_run << "cur_dir = os.getcwd() " << std::endl;
		os_run << "PLATO_root  = \"" << plato_root_dir << "\"" << std::endl;
		os_run << "PLATO_path = os.path.join(PLATO_root,\"PLATO\")" << std::endl << std::endl;

		os_run << "os.chdir(PLATO_root)" << std::endl << std::endl;

		os_run << "# Load all the files related to the data analysis" << std::endl;
		os_run << "execfile(os.path.join(PLATO_path,\"load_PLATO_files.py\"))" << std::endl << std::endl;

		os_run << "# Load all the files to manipulate filenames / file content / ...; create files / directories etc" << std::endl;
		os_run << "execfile(os.path.join(PLATO_path,\"load_files_for_XML.py\"))" << std::endl << std::endl;

		os_run << "os.chdir(cur_dir) " << std::endl << std::endl;

		os_run << "xml2PLATO(\"" << plato_input_fname << "\")" << std::endl;
	}
	catch( std::exception& ex )
	{
		PrintLog("Error in CollectCrdAnalMod::SavePlatoInputFiles() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int CollectCrdAnalMod::RunPlato()
{
	PrintLog(" Run PLATO program \n");
	int ires = pApp->ExecuteScriptFromFile(plato_run_fname.c_str());
	return ires;
}

int CollectCrdAnalMod::CalcTimeProj()
{
	char msg[256];
	char* buf = NULL;
	try
	{
		AtomGroup* p_at_grp = GetActiveAtomGroup();
		if( p_at_grp == NULL ) throw std::runtime_error("No Active Atom Group");
		int na = p_at_grp->GetNAtoms();
		if( na == 0 ) throw std::runtime_error(" Empty Active Atom Group ");

		int nv = eigen_vecs.size();
		if( nv == 0 ) throw std::runtime_error(" No Eigen Vectors set ");
		int iv;
		int ncrd;
		for( iv = 0; iv < nv; iv++)
		{
			ncrd = eigen_vecs[iv].size();
			if( ncrd != 3*na )
			{
				sprintf(msg,"Num coords of Eigen vector %d is %d  not equal to Num coords of active atom group %d",
					        iv, ncrd, 3*na);
				throw std::runtime_error(msg);
			}
		}

		time_projections.clear();
		time_projections.resize(nv);

		int bsize = 60*na;
		buf = (char*) malloc( bsize );
		
		ifstream is(md_traj_fname_plato.c_str());
		if( is.fail()) throw std::runtime_error("Error to open file " + md_traj_fname_plato);

		is.getline(buf, bsize);
		for(;;)
		{
			is.getline(buf, bsize);
			if( is.fail() ) break;
			std::vector<std::string> str_arr;
			std::string line = buf;
			boost::split(str_arr,line,boost::is_any_of(";"),boost::token_compress_on);
			ncrd = str_arr.size();
			if( ncrd != 3*na )
			{
				sprintf(msg,"Num coords in traj file line %d is not equal to Num coords of active atom group %d",
					        ncrd,3*na);
				throw std::runtime_error(msg);
			}
			int j;
			HaVec_double crd(3*na);
			for( j=0; j < 3*na; j++)
			{
				crd[j] = atof(str_arr[j].c_str()); 
			}
			for(iv = 0; iv < nv; iv++)
			{
				double dp = dot_product(crd,eigen_vecs[iv]);
				time_projections[iv].push_back(dp);
			}
		}
	}
	catch( std::exception& ex)
	{
		if( buf != NULL ) free (buf);
		PrintLog("Error in CollectCrdAnalMod::CalcTimeProj() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int CollectCrdAnalMod::SaveTimeProjFile()
{
	char buf[256];
	try
	{
		std::ofstream os(time_proj_fname.c_str());
		if( os.fail() ) throw std::runtime_error("Error to open file " + time_proj_fname );
        
		int nv  = time_projections.size();
		if( nv == 0 ) throw std::runtime_error("No time projections " );
		int npt = time_projections[0].size();
		int i,iv;
		for( iv = 0; iv < nv; iv++ )
		{
			sprintf(buf," %6d       ",(iv+1));
			os << buf;
		}
		os << std::endl;

		for( i = 0; i < npt; i++ )
		{
			for( iv = 0; iv < nv; iv++)
			{
				double val = 99999.99;
				if( i < time_projections[iv].size() ) val = time_projections[iv][i];
				sprintf(buf," %12.6f ",val);
				os << buf;
			}
			os << std::endl;
		}
	}
	catch( std::exception& ex )
	{
		PrintLog("Error in CollectCrdAnalMod::SaveTimeProjFile() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;	
}

int CollectCrdAnalMod::ShiftAlongEigenVec( int idx_vec, double shift_val )
{
//	PrintLog(" CollectCrdAnalMod::ShiftAlongEigenVec() pt 1  idx_vec = %d shift_val=%9.3f \n", idx_vec, shift_val );
	try
	{
		if( idx_vec < 0 || idx_vec >= eigen_vecs.size() ) throw std::runtime_error("invalid eigenvector index ");
		AtomGroup* pgrp = GetActiveAtomGroup();
		if( pgrp == NULL ) throw std::runtime_error("Active Atom Group is not set ");
		if( pgrp->GetNAtoms()*3 != eigen_vecs[idx_vec].size() ) throw std::runtime_error("Dimension of the eigen vector is invalid ");

		AtomIteratorAtomGroup aitr(pgrp);
		HaAtom* aptr;

		int i = 0;
		for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
		{
			aptr->SetX( aptr->GetX() + shift_val* eigen_vecs[idx_vec][3*i] ) ;
			aptr->SetY( aptr->GetY() + shift_val* eigen_vecs[idx_vec][3*i+1] ) ;
			aptr->SetZ( aptr->GetZ() + shift_val* eigen_vecs[idx_vec][3*i+2] ) ;
			i++;
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in CollectCrdAnalMod::ShiftAlongEigenVec() \n");
		PrintLog(" %s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

