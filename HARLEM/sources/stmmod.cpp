/*! \file stmmod.cpp

    Classes to model STM current 

    \author Igor Kurnikov  
    \author Jianping, Duke University 
 
    \date 2001-2002
*/

#define STMMOD_CPP

#include <mpi.h>

#include "stmmod.h"
#include "haqchem.h"
#include "hamolset.h"
#include "hamultipole.h"
#include "math.h"

StmMod::StmMod(HaMolSet* new_phost_mset):
HaCompMod(COMP_MOD_STM,new_phost_mset)
{
	SetStdParams();
}

StmMod::~StmMod()
{
	
}

void
StmMod::SetStdParams()
{
    num_mol_ao = 56;
	tmat_el = 0.0;
	wave_ve_1 = 0.0;
	wave_ve_2 = 0.0;
	wave_ve_3 = 0.0;
}

int
StmMod::CalcTMatr1()
{
	HaMolSet* pmset = GetMolSet();

	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
	if(ptr_qc_mod == NULL)
	{
		ErrorInMod("StmMod::CalcTMatr1()", "QChem module is not initialized");
		return FALSE;
	}

    int nb = ptr_qc_mod->GetNBfunc(); // Get Number of basis functions

	cout << "Print 1 " << endl;

//	ovlp.Print_info(cout,1);
	cout << endl;

	cout << " Print 2 " << endl;

	//HaOperR r1(*ptr_qc_mod);   // Initial Electric Dipole  operator
	//r1.Print_info(cout,1);	
//cout << endl;

    //HaOperGrad rg1(*ptr_qc_mod);  // Initiate gradient operator
	//rg1.Print_info(cout,1);	
	//cout << endl;

	//HaOperRDelt rd1(*ptr_qc_mod);  // Initiate R X Delt (magnetic dipole) operator
//	cout << " r X Grad for non-london orbitals";
	//rd1.Print_info(cout,1);

//	rd1.RecalcLondon();
//	cout << " r X Grad for london orbitals";
//	rd1.Print_info(cout,1);	
//	cout << endl;

	HaMat_double ss,hh,sp; // Fill Overlap Matrix
	
	ss = ptr_qc_mod->GetOvlpMat();
	int i,j,nev,m,n,a,b;
	hh.newsize(nb,nb);
	sp.newsize(nb,nb);
	//for(i=1; i<=nb;i++)
	//{
		//for(j=1;j<=nb;j++)
		//{
			//hh(j,i)=1.75*ss(j,i);
			//cout<< hh(j,i)<<endl;
	//	}
	//cout<<"\n";	
	//}

    
	//HaMat_double rx,ry,rz;
	//HaMat_double rmx,rmy,rmz;
	

	//r1.FillMat(rx,1); // Fill matricies rx, ry, rz with matricies of electric dipole operator  
//	r1.FillMat(ry,2);
	//r1.FillMat(rz,3);
//
//	rd1.FillMat(rmx,1); // Fill matricies rx, ry, rz with matricies of magnetic dipole operator  
//	rd1.FillMat(rmy,2);
//	rd1.FillMat(rmz,3);


	//FILE* fsmat = fopen("overlap_mat.txt","w"); // Save overlap matrix into overlap_mat.tat file 
	//int i, j;

	//if(fsmat != NULL)
	//{
	//	fprintf(fsmat,"%5d%5d\n",nb,nb);
	//	for(i = 1; i <= nb; i++)
	//	{
		//	const HaAtom* aptr = ptr_qc_mod->GetAtomOfAO(i);
		//	int elem= aptr->GetElemNo();

		//	fprintf(fsmat,"%4d ",elem);

		//	for(j = 1; j <= nb; j++)
		//	{
		//		//hh(j,i)=1.75*ss(j,i);
				//fprintf(fsmat,"%16.9f ",ss(j,i));
		//		//fprintf(fsmat,"%16.9f",hh(j,i));
		//	}
		//			fprintf(fsmat,"\n");
		//}
	//}
	//fclose(fsmat);

    HaMat_double sd,sh,shs,eigv,rev;
	HaVec_double eig;
	//HaVec_double ip(nb);

	//ip(1) = 10.0;
//	ip(6) =  12.0;
	

	//HaMat_double::mat_sdiag(ss,sd,eig);
	sp = ss;
	sp.SqRoot(-1);

//	cout << " Eigen Values of SS matrix" << endl;
//	for( i = 1; i <= nb; i++)
//	{
	   
//		const HaAtom* aptr = ptr_qc_mod->GetAtomOfAO(i);
//		int elem= aptr->GetElemNo();
		
//		if( elem == 6 )
//		{
//        hh(i,j)=-21.4;
//		 hh(i+1,j+1)=-11.4;
//		 hh(i+2,j+2)=-11.4;
//		 hh(i+3,j+3)=-11.4;
//		 i+=3;
//		}
//	else if( elem == 1)
//		{
 //         hh(i,j)=-13.6;
	//	}
//		else if( elem == 7)
//		{
                            

//		}
	
//		cout << i << "  " << eig(i) << endl;
	   
//	}
	cout << " Huckel Matrix" << endl;
	//no=0;
	nev=0;
	for( i = 1; i <= nb; i++)
	{
		const HaAtom* aptr = ptr_qc_mod->GetAtomOfAO(i);
		int elem= aptr->GetElemNo();
		
		if( elem == 6 )
		{
         hh(i,i)=-21.4;
		 hh(i+1,i+1)=-11.4;
		 hh(i+2,i+2)=-11.4;
		 hh(i+3,i+3)=-11.4;
		 i+=3;
		 //no+=4;
		 //nev+=4;
		}
		
	else if( elem == 1)
		{
          hh(i,i)=-13.6;
		  //no+=1;
		  //nev+=1;
		}
		else if( elem == 7)
		{
         hh(i,i)=-26.0;
		 hh(i+1,i+1)=-13.4;
		 hh(i+2,i+2)=-13.4;
		 hh(i+3,i+3)=-13.4;
		 i+=3;
		 //no+=4;
		 //nev+=5;
		}
		
	else if( elem==8)
	{
         hh(i,i)=-32.3;
		 hh(i+1,i+1)=-14.8;
		 hh(i+2,i+2)=-14.8;
		 hh(i+3,i+3)=-14.8;
		 i+=3;
		 //no+=4;
		 //nev+=6;
	}
	else if( elem==14)
	{
         hh(i,i)=-17.3;
		 hh(i+1,i+1)=-9.2;
		 hh(i+2,i+2)=-9.2;
		 hh(i+3,i+3)=-9.2;
		 i+=3;
		 //no+=4;
		 //nev+=4;
	}
		else if( elem==15)
	{
         hh(i,i)=-18.6;
		 hh(i+1,i+1)=-14.0;
		 hh(i+2,i+2)=-14.0;
		 hh(i+3,i+3)=-14.0;
		 i+=3;
		 //no+=4;
		 //nev+=5;
	}
	else if( elem==16)
	{
         hh(i,i)=-20.0;
		 hh(i+1,i+1)=-13.3;
		 hh(i+2,i+2)=-13.3;
		 hh(i+3,i+3)=-13.3;
		 i+=3;
		 //no+=4;
		 nev+=6;
	}
	else if( elem==79)
	{
         hh(i,i)=-10.92;
		 //hh(i+1,i+1)=-5.55;
		 //hh(i+2,i+2)=-5.55;
		 //hh(i+3,i+3)=-5.55;
		 //hh(i+4,i+4)=-15.07;
         //hh(i+5,i+5)=-15.07;
         //hh(i+6,i+6)=-15.07;
         //hh(i+7,i+7)=-15.07;
         //hh(i+8,i+8)=-15.07;
		 //i+=8;
		 //no+=4;
		// nev+=1;
	}
	else if( elem==29)
	{
         hh(i,i)=-11.4;
		 hh(i+1,i+1)=-6.06;
		 hh(i+2,i+2)=-6.06;
		 hh(i+3,i+3)=-6.06;
		 hh(i+4,i+4)=-14.0;
         hh(i+5,i+5)=-14.0;
         hh(i+6,i+6)=-14.0;
         hh(i+7,i+7)=-14.0;
         hh(i+8,i+8)=-14.0;
		 i+=8;
		 //no+=4;
		// nev+=11;
	}
		//cout << i << "  " << eig(i) << endl;
	   
	}
//cout<<"Orbital number: "<<no<<endl;
//cout <<"Number of valent electron: "<<nev<<endl;
    //hh(1,1)=-4.880;
	hh(num_mol_ao+1,num_mol_ao+1)=-4.86;
	for (i=1;i<=nb;i++)
	{
		for (j=1;j<=nb;j++)
		{if(i!=j)
				hh(i,j)=1.75*ss(i,j)*(hh(i,i)+hh(j,j))/2.0;
		}
	}

    

/*	char buf[256];

	for(i=1;i<=nb;i++)
	{
		for (j=1;j<=nb;j++)
		{
			sprintf(buf, "%7.4f ",hh(i,j));
          //sprintf(buf, "%7.4f ",sp(i,j));
			cout << buf ;
		}
		cout << endl;
	}*/
//	cout << hh(1,2) <<endl;
     matmult(sh,sp,hh);
     matmult(shs,sh,sp);
     HaMat_double::mat_sdiag(shs,eigv,eig);
	 matmult(rev,sp,eigv);
//cout<< " Eigenvectors"<<endl;
//char buf[256];
for(i=1;i<=nb;i++)
	{
		/*for (j=1;j<=nb;j++)
		{
			//sprintf(buf, "%7.4f ",shs(i,j));
          sprintf(buf, "%7.4f ",rev(i,j));
		  cout << buf ;
		}
		cout << endl;*/
	cout <<i<<" "<<eig(i)<<endl;

}

/*cout << hh(2,2) << endl;
cout << hh(5,5) << endl;
cout << ss(1,5) << endl;
cout << ss(2,5) << endl;
cout << ss(3,5) << endl;
cout << ss(4,5) << endl;
cout << ss(2,6) << endl;
cout << ss(3,6) << endl;
cout << ss(4,6) << endl;
cout << ss(3,7) << endl;
cout << ss(4,7) << endl;
cout << ss(4,8) << endl; */
// Density of states
/*int counter;
double estart,emax,emin;
for(s=1;s<=54;s++)                              
	{
		estart=-16.0+1.0*(s-1);
		emax=estart+0.5;
		emin=estart-0.5;
		counter=0;
		for(i=1;i<=nb;i++)
			{
				if(eig(i)>=emin && eig(i)<=emax)
						counter++;
			}
		cout <<estart << " " <<counter<<endl;
	}*/			
//double sum,sum1,sum2,kx,ky,kz,la;
//sum1=0;
//sum2=0;
//la=2.884;
//kx=3.1415926/(2*la);
//ky=3.1415926/(2*la);
//kz=3.1415926/(2*la);
//l=2;
//for(l=1;l<=3;l++)
//{
/*for(k=0;k<=100;k++)
	{
	sum1=0;
	sum2=0;
	double e0,vc,e;
	e=-5.30;
	e0=-10.92;
	vc=-6.666;
		for(m=1;m<=9;m++)
			{
				for(n=1;n<=9;n++)
					{
					  sum1=sum1+cos(m*k*3.1415926/100+n*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)))*hh(82,9*m-n+1);
					  sum2=sum2+sin(m*k*3.1415926/100+n*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)))*hh(82,9*m-n+1);
					//sum1=sum1+cos(kx*m*la+ky*n*la+kz*l*la)*hh(301,100*(l-1)+10*m-n+1);
					//sum2=sum2+sin(kx*m*la+ky*n*la+kz*l*la)*hh(301,100*(l-1)+10*m-n+1);
					//sum1=sum1+cos(kx*m*la+kx*n*la/2+ky*n*1.73205*la/2)*hh(101,10*m-n+1);
					//sum2=sum2+sin(kx*m*la+kx*n*la/2+ky*n*1.73205*la/2)*hh(101,10*m-n+1);
					//sum1=sum1+cos(k*m*3.1415926/100+k*n*3.1415926/200+n*acos(((e-e0)/vc-2*cos(k*3.1415926/100))/(4*cos(k*3.1415926/200))))*hh(101,10*m-n+1);
					//sum2=sum2+sin(k*m*3.1415926/100+k*n*3.1415926/200+n*acos(((e-e0)/vc-2*cos(k*3.1415926/100))/(4*cos(k*3.1415926/200))))*hh(101,10*m-n+1);
					}
			}
	
	//}

	sum=sum1*sum1+sum2*sum2;
	cout << sum <<endl;
	}*/
//cout << sum2 <<endl;					
/*cout<< "Transimission Matrix"<<endl;
HaMat_double sg,sgs,esb,greenb;
double ga1,ga2,ga3,ga4;
mat_scale(esb,ss,-4.86);
mat_diff(greenb,esb,hh);
HaMat_double::mat_inverse(greenb);*/
//matmult(sg,ss,greenb);
//matmult(sgs,sg,ss);
//ga1=greenb(49,54)*greenb(49,54);
//tmat_el=ga1;
//ga2=greenb(83,87)*greenb(83,87);
//wave_ve_1=ga2;
//ga3=greenb(84,88)*greenb(84,88);
//wave_ve_2=ga3;
//ga4=greenb(25,32)*greenb(25,32);
//wave_ve_3=ga4;
//cout << sgs(1,189) << endl;
//cout << sgs(2,190) << endl;
//cout << sgs(3,191) << endl;
//cout << sgs(4,192) << endl;

//cout << sgs(2,nb-6)<< endl;
//cout << sgs(3,nb-5) << endl;
//cout << sgs(4,nb-4) << endl;
//cout << sgs(5,nb-3) << endl;
//cout << " ---------------------------------------------"<<endl;

//HaMat_double es[50],green[50],sg[50],sgs[50];
/*HaMat_double ng[1];
ng[0].newsize(nb,nb);
ng[1].newsize(nb,nb);
ng[2].newsize(nb,nb);
ng[3].newsize(nb,nb);
ng[4].newsize(nb,nb);*/
/*for(k=0;k<=0;k++)
{
	mat_scale(es[k],ss,(eig(nev/2)+(k+1)*(eig(nev/2+1)-eig(nev/2))/6));*/
	//mat_scale(es[k],ss,(-10.0+(k+1)*9.0/51));
	/*mat_diff(green[k],es[k],hh);
	HaMat_double::mat_inverse(green[k]);
	matmult(sg[k],ss,green[k]);
	matmult(sgs[k],sg[k],ss);*/
//	double sum;
/*	double sum,sum1,sum2;
	int x;
	char buf[256];
	for(l=1;l<=nb;l++)
	{
			for (m=1;m<=nb;m++)
		{
		    sum=0;
			for(j=1;j<=nb;j++)
			{
			//	sum=sum+rev(l,j)*rev(j,m)/(eig(nev/2)+(k+1)*(eig(nev/2+1)-eig(nev/2))/6-eig(j));
				sum1=0;
				sum2=0;
				for (x=1; x<=nb;x++)
				{
					sum1=sum1+rev(x,j)*ss(l,x);
					sum2=sum2+rev(x,j)*ss(x,m);
				}
				//sum=sum+sum1*sum2/(eig(nev/2)+(k+1)*(eig(nev/2+1)-eig(nev/2))/6-eig(j));
				sum=sum+sum1*sum2/(-5.30-eig(j));
			}
			ng[0](l,m)=sum;
            sprintf(buf, "%16.9f ",ng[0](l,m));
			cout << buf;
		}
		cout<<endl;

		}*/


	/*char buf[256];

	for(i=1;i<=nb;i++)
	{
		for (j=1;j<=nb;j++)
		{
			//sprintf(buf, "%16.9f ",green[k](i,j));
           sprintf(buf, "%16.9f ",sgs[k](i,j));
			cout << buf ;
		}
		cout <<endl;
	
	}
	cout<<"----------------------"<<endl;*/
	//cout << (1-((-3.0+(k+1)*6.0/3)/(-2*0.05882*27.2))*((-3.0+(k+1)*6.0/3)/(-2*0.05882*27.2)))*sgs[k](25,29)*sgs[k](25,29) << endl;
	//cout << sgs[k](33,37) << endl;
//	cout << sgs[k](17,174)/sgs[k](17,173)+sgs[k](17,173)/sgs[k](17,174)<<endl;
//	cout<<sgs[k](17,174)<<sgs[k](17,175)<<endl;
	//cout<<"------------------------------------------------------------------"<<endl;
//}
//for (k=0;((k*6+9)<=nb);k++)
//cout << sgs(2,9+k*6)<<endl;
	/*sum=0;
	for(s=1;s<=25;s++)
	{
		for(a=27;a<=nb;a++)
			{
				for(b=27;b<=nb;b++)
					{
						sum=sum+(hh(s,a)-(-10.0+(k+1)*9.0/50)*ss(s,a))*sgs[k](a,b)*(hh(b,26)-(-10.0+(k+1)*9.0/3)*ss(b,26));
					}
			}
	}

	cout<<sum*sum << endl;*/
//}
//HaMat_double es[50],green[50],sg[50],sgs[50];
//double tot,tot1,tot2,e0,vc,e;
//e=-5.30;
//e0=-0.29456*27.2;
//vc=-0.05882*27.2;
//for(k=0;k<=99;k++)
//{
	//HaMat_double es,green,sg,sgs;
//	mat_scale(es,ss,(e0+2*vc*(cos(kix*lai)+sin(kiy*lai));
//	mat_scale(es,ss,e);
//	mat_scale(es,ss,(-10.0+(k+1)*9.0/51));
	/*mat_scale(es,ss,e0);
	mat_diff(green,es,hh);
	HaMat_double::mat_inverse(green);
	matmult(sg,ss,green);
	matmult(sgs,sg,ss);
	tot1=0;
	tot2=0;*/
//	int mod;
	/*l=5;
	for(s=1;s<=l*l;s++)
		{
			for(a=l*l+2;a<=nb;a++)
				{
					for(b=l*l+2;b<=nb;b++)
						{
							//tot1=tot1+(hh(s,a)-e*ss(s,a))*sgs(a,b)*(hh(b,101)-e*ss(b,101))*cos((s+10-(s%10))*k*3.1415926/(10*100)+(11-(s%10))*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)));
                            //tot2=tot2+(hh(s,a)-e*ss(s,a))*sgs(a,b)*(hh(b,101)-e*ss(b,101))*sin((s+10-(s%10))*k*3.1415926/(10*100)+(11-(s%10))*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)));
							tot1=tot1+(hh(s,a)-e*ss(s,a))*sgs(a,b)*(hh(b,l*l+1)-e*ss(b,l*l+1))*cos((s+l-(s%l))*3.1415926/(l*2)+(l+1-(s%l))*3.1415926/2);
                            tot2=tot2+(hh(s,a)-e*ss(s,a))*sgs(a,b)*(hh(b,l*l+1)-e*ss(b,l*l+1))*sin((s+l-(s%l))*3.1415926/(l*2)+(l+1-(s%l))*3.1415926/2);
						}
				}
		}
	tot=tot1*tot1+tot2*tot2;
	cout << tot<< endl;*/
//	cout << s << endl;
//}
/*	double tot,tot1,tot2,e0,vc,e;
	int q;
	e=-5.30;
	e0=-10.92;
	vc=-6.666;
	for(k=30;k<=30;k++)
		{
			HaMat_double es,green,sg,sgs;*/
//	mat_scale(es,ss,(e0+2*vc*(cos(kix*lai)+sin(kiy*lai));
		//	mat_scale(es,ss,e);
//	mat_scale(es,ss,(-10.0+(k+1)*9.0/51));
//	mat_scale(es,ss,e0);
/*	mat_diff(green,es,hh);
	HaMat_double::mat_inverse(green);
	matmult(sg,ss,green);
	matmult(sgs,sg,ss);
	tot1=0;
	tot2=0;
	q=10;*/
//	for(l=1;l<=4;l++)
//		{
		/*	for(m=1;m<=q;m++)
				{
					for(n=1;n<=q;n++)
						{
							for(a=(q*q+2);a<=nb;a++)
								{
									for(b=(q*q+2);b<=nb;b++)
										{
											tot1=tot1+(hh(q*m-n+1,a)-e*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-e*ss(b,q*q+1))*cos(m*k*3.1415926/100+n*k*3.1415926/200+n*acos(((e-e0)/vc-2*cos(k*3.1415926/100))/(4*cos(k*3.1415926/200))));
											tot2=tot2+(hh(q*m-n+1,a)-e*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-e*ss(b,q*q+1))*sin(m*k*3.1415926/100+n*k*3.1415926/200+n*acos(((e-e0)/vc-2*cos(k*3.1415926/100))/(4*cos(k*3.1415926/200))));
										  //tot1=tot1+(hh(10*m-n+1,a)-e*ss(10*m-n+1,a))*sgs(a,b)*(hh(b,101)-e*ss(b,101))*cos(m*k*3.1415926/100+n*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)));
										  //tot2=tot2+(hh(10*m-n+1,a)-e*ss(10*m-n+1,a))*sgs(a,b)*(hh(b,101)-e*ss(b,101))*sin(m*k*3.1415926/100+n*acos((e-e0)/(2*vc)-cos(k*3.1415926/100)));
										  //tot1=tot1+(hh(100*(l-1)+10*m-n+1,a)-e*ss(100*(l-1)+10*m-n+1,a))*sgs(a,b)*(hh(b,401)-e*ss(b,401))*cos(m*3.1415926/2+n*3.1415926/2+l*3.1415926/2);
										  //tot2=tot2+(hh(100*(l-1)+10*m-n+1,a)-e*ss(100*(l-1)+10*m-n+1,a))*sgs(a,b)*(hh(b,401)-e*ss(b,401))*sin(m*3.1415926/2+n*3.1415926/2+l*3.1415926/2);
										}
								}
						}
				}*/
//		}
	//tot=tot1*tot1+tot2*tot2;
//	cout << tot<< endl;
//	cout << s << endl;
	//	}
//cout <<sgs(90,107)<<endl;
/*cout << " print 5" << endl;
	double tot,tot1,tot2,estart,emax,emin;
	double ener(int i,int j,int k);
	int q,counter,v;
	for(v=1;v<=22000;v++)
		{
			estart=-4.330+0.001*v;
			emax=estart+0.0005;
			emin=estart-0.0005;
			counter=0;
			tot=0;
			for(i=1;i<=11;i++)
			{
			for(j=1;j<=11;j++)
			{
			for(k=1;k<=10;k++)
			{
			tot1=0;
	        tot2=0;
			if(ener(i,j,k)>=emin && ener(i,j,k)<emax)
			{
			HaMat_double es,green,sg,sgs;
	        mat_scale(es,ss,ener(i,j,k));
        	mat_diff(green,es,hh);
	        HaMat_double::mat_inverse(green);
	        matmult(sg,ss,green);
	        matmult(sgs,sg,ss);
	        q=11;
			for(m=1;m<=q;m++)
			{
			for(n=1;n<=q;n++)
			{
			for(a=(q*q+2);a<=nb;a++)
			{
			for(b=(q*q+2);b<=nb;b++)
			{
			tot1=tot1+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*cos(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			tot2=tot2+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*sin(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			//tot1=tot1+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*cos(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			//tot2=tot2+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*sin(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			}
			}
			}
			}
			counter++;
			}
			tot=tot+tot1*tot1+tot2*tot2;
			}
			}
			}
			if(counter>0)
     cout << tot <<"                "  <<v << "             " << counter << endl;
     }
}

double ener(int a,int b,int c)
{
int qe=2*11;
return (-10.920-4*6.666*(cos(2*3.1415926*b/qe-2*3.1415926*a/qe+2*3.1415926*c/20)*cos(-2*3.1415926*b/qe+2*3.1415926*a/qe+2*3.1415926*c/20)+
cos(-2*3.1415926*b/qe+2*3.1415926*a/qe+2*3.1415926*c/20)*cos(2*3.1415926*b/qe+2*3.1415926*a/qe-2*3.1415926*c/20)+
cos(2*3.1415926*b/qe+2*3.1415926*a/qe-2*3.1415926*c/20)*cos(2*3.1415926*b/qe-2*3.1415926*a/qe+2*3.1415926*c/20)));	
}*/

/*cout << " print 5" << endl;
	double tot,tot1,tot2,estart,emax,emin;
	double ener(int i,int j,int k);
	int q,counter,v;
	//tot=0;
	for(v=1;v<=30;v++)
		{
			estart=-4.883+0.05*v;
			//emin=-10.92+0.1*(v-1);
			counter=0;
			tot=0;
			for(i=1;i<=19;i++)
			{
			for(j=1;j<=19;j++)
			{
			for(k=1;k<=10;k++)
			{
			tot1=0;
	        tot2=0;
			if(ener(i,j,k)>-4.880 && ener(i,j,k)<=estart)
			{
			HaMat_double es,green,sg,sgs;
	        mat_scale(es,ss,ener(i,j,k));
        	mat_diff(green,es,hh);
	        HaMat_double::mat_inverse(green);
	        matmult(sg,ss,green);
	        matmult(sgs,sg,ss);
	        q=19;
			for(m=1;m<=q;m++)
			{
			for(n=1;n<=q;n++)
			{
			for(a=(q*q+2);a<=nb;a++)
			{
			for(b=(q*q+2);b<=nb;b++)
			{
			tot1=tot1+(hh(q*(n-1)+m,a)-ener(i,j,k)*ss(q*(n-1)+m,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*cos(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			tot2=tot2+(hh(q*(n-1)+m,a)-ener(i,j,k)*ss(q*(n-1)+m,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*sin(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			//tot1=tot1+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*cos(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			//tot2=tot2+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*sin(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			}
			}
			}
			}
			counter++;
			}
			tot=tot+tot1*tot1+tot2*tot2;
			}
			}
			}
			if(counter>0)
     cout << tot <<"                "  <<v << "             " << counter << endl;
     }
}

double ener(int a,int b,int c)
{
int qe=2*19;
return (-10.920-4*6.666*(cos(2*3.1415926*b/qe-2*3.1415926*a/qe+2*3.1415926*c/20)*cos(-2*3.1415926*b/qe+2*3.1415926*a/qe+2*3.1415926*c/20)+
cos(-2*3.1415926*b/qe+2*3.1415926*a/qe+2*3.1415926*c/20)*cos(2*3.1415926*b/qe+2*3.1415926*a/qe-2*3.1415926*c/20)+
cos(2*3.1415926*b/qe+2*3.1415926*a/qe-2*3.1415926*c/20)*cos(2*3.1415926*b/qe-2*3.1415926*a/qe+2*3.1415926*c/20)));	
}*/

//cout << hh(1,6) << "   "  << ss(1,6) <<endl;


 /*HaMat_double ng,nss;
     ng.newsize(248,248);
	 nss.newsize(248,248);
	 for(i=1;i<=248;i++)
	 {
	 for(j=1;j<=248;j++)
	 {
	 ng(i,j)=hh(i+1,j+1);
	 nss(i,j)=ss(i+1,j+1);
	 }
	 }

cout << " print 5" << endl;
	double tot,tot1,tot2,estart,emax,emin;
	double ener(int i,int j,int k);
	int q,counter,v;
	//tot=0;
	//counter=0;
	for(v=1;v<=1;v++)
		{
			//estart=-4.880+0.05*v;
			estart=-5.047+0.001*v;
			emin=estart-0.0005;
			emax=estart+0.0005;
			counter=0;
			tot=0;
			//tot1=0;
			//tot2=0;
			for(i=1;i<=19;i++)
			{
			for(j=1;j<=19;j++)
			{
			for(k=1;k<=10;k++)
			{
			tot1=0;
	        tot2=0;
			if(ener(i,j,k)>emin && ener(i,j,k)<=emax)
			{
			HaMat_double es,green,sg,sgs;
	        mat_scale(es,nss,ener(i,j,k));
        	mat_diff(green,ng,es);
	        HaMat_double::mat_inverse(green);
	       // matmult(sg,ss,green);
	       // matmult(sgs,sg,ss);
	        q=19;
			for(m=1;m<=q;m++)
			{
			for(n=1;n<=q;n++)
			{
			for(a=2;a<=nb-q*q;a++)
			{
			for(b=2;b<=nb-q*q;b++)
			{
			//tot1=tot1+(hh(nb-q*q+q*(n-1)+m,a)-ener(i,j,k)*ss(nb-q*q+q*(n-1)+m,a))*green(a,b)*(hh(b,57)-ener(i,j,k)*ss(b,57))*cos(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			//tot2=tot2+(hh(nb-q*q+q*(n-1)+m,a)-ener(i,j,k)*ss(nb-q*q+q*(n-1)+m,a))*green(a,b)*(hh(b,57)-ener(i,j,k)*ss(b,57))*sin(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			tot1=tot1+(hh(nb-q*q+q*(n-1)+m,a)-ener(i,j,k)*ss(nb-q*q+q*(n-1)+m,a))*green(a-1,b-1)*(hh(b,1)-ener(i,j,k)*ss(b,1))*cos(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			tot2=tot2+(hh(nb-q*q+q*(n-1)+m,a)-ener(i,j,k)*ss(nb-q*q+q*(n-1)+m,a))*green(a-1,b-1)*(hh(b,1)-ener(i,j,k)*ss(b,1))*sin(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			//tot1=tot1+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*cos(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			//tot2=tot2+(hh(q*m-n+1,a)-ener(i,j,k)*ss(q*m-n+1,a))*sgs(a,b)*(hh(b,q*q+1)-ener(i,j,k)*ss(b,q*q+1))*sin(2*3.1415926*i*m/21+2*3.1415926*j*n/21);
			}
			}
			}
			}
			tot=tot+tot1*tot1+tot2*tot2;
			counter++;
			}
		//	tot=tot+tot1*tot1+tot2*tot2;
			}
			}
			}
			//tot=tot+tot1*tot1+tot2*tot2;
			if(counter>0)
     cout << tot <<"                "  <<v << "             " << counter << endl;
     }*/
	 /*HaMat_double ng,nss;
     ng.newsize(56,56);
	 nss.newsize(56,56);
	 for(i=1;i<=56;i++)
	 {
	 for(j=1;j<=56;j++)
	 {
	 ng(i,j)=hh(i+1,j+1);
	 nss(i,j)=ss(i+1,j+1);
	 }
	 }*/
	        /*HaMat_double es,green,sg,sgs;
	        mat_scale(es,ss,-4.86);
        	mat_diff(green,hh,es);
	        HaMat_double::mat_inverse(green);
			//matmult(sg,ss,green);
	        //matmult(sgs,sg,ss);
			cout << green(1,37)*green(1,37) << "  " << green(1,37)<< endl;
}*/

/*double ener(int a,int b,int c)
{
double kxa,kya,kza;
int qe=2*19;
kxa=2*2*3.1415926*a/qe;
kya=2*(2*3.1415926*b*2*1.7320508/(3*qe)-2*3.1415926*a*1.7320508/(3*qe));
kza=2*(2*3.1415926*a*2.44948974/(6*qe)+2*3.1415926*b*2.44948974/(6*qe)-2*3.1415926*c*2.44948974/(4*10));
return (-10.920-2*6.666*(cos(kxa)+cos(kxa/2+1.7320508*kya/2)+cos(kxa/2-1.7320508*kya/2)+
cos(kxa/2+1.7320508*kya/6+2.44948974*kza/3)+cos(1.7320508*kya/3-2.44948974*kza/3)+
cos(kxa/2-1.7320508*kya/6-2.44948974*kza/3)));
}*/


/*cout << " print 5" << endl;
	double tot,tot1,tot2,estart,emax,emin;
	double ener(int i,int j,int k);
	int q,counter,v;
		for(v=1;v<=1;v++)
		{
			estart=-5.047+0.001*v;
			emin=estart-0.0005;
			emax=estart+0.0005;
			counter=0;
			tot=0;
			for(i=1;i<=19;i++)
			{
			for(j=1;j<=19;j++)
			{
			for(k=1;k<=10;k++)
			{
			tot1=0;
	        tot2=0;
			if(ener(i,j,k)>emin && ener(i,j,k)<=emax)
			{
	        q=19;
			for(m=1;m<=q;m++)
			{
			for(n=1;n<=q;n++)
			{
			tot1=tot1+hh(1+q*(n-1)+m,1)*cos(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			tot2=tot2+hh(1+q*(n-1)+m,1)*sin(2*3.1415926*i*m/q+2*3.1415926*j*n/q+2*3.1415926*k*1/10);
			}
			}
			tot=tot+tot1*tot1+tot2*tot2;
			counter++;
			}
			}
			}
			}
			if(counter>0)
     cout << tot <<"                "  <<v << "             " << counter << endl;
     }
}

double ener(int a,int b,int c)
{
double kxa,kya,kza;
int qe=2*19;
kxa=2*2*3.1415926*a/qe;
kya=2*(2*3.1415926*b*2*1.7320508/(3*qe)-2*3.1415926*a*1.7320508/(3*qe));
kza=2*(2*3.1415926*a*2.44948974/(6*qe)+2*3.1415926*b*2.44948974/(6*qe)-2*3.1415926*c*2.44948974/(4*10));
return (-10.920-2*6.666*(cos(kxa)+cos(kxa/2+1.7320508*kya/2)+cos(kxa/2-1.7320508*kya/2)+
cos(kxa/2+1.7320508*kya/6+2.44948974*kza/3)+cos(1.7320508*kya/3-2.44948974*kza/3)+
cos(kxa/2-1.7320508*kya/6-2.44948974*kza/3)));
}*/

     int z;
	 z= num_mol_ao;
	 HaMat_double ng,nss;
     ng.newsize(z,z);
	 nss.newsize(z,z);
	 for(i=1;i<=z;i++)
	 {
	 for(j=1;j<=z;j++)
	 {
	 ng(i,j)=hh(i,j);
	 nss(i,j)=ss(i,j);
	 }
	 }
	double tot,tot1,tot2,enr,k1,k2,k3,pi;
	int q;
			enr=-4.86;
			pi=3.1415926;
		//	k1=-0.494;
		//	k2=0.409;
		//	k3=0.348;
			k1=wave_ve_1;
			k2=wave_ve_2;
			k3=wave_ve_3;
			tot1=0;
	        tot2=0;
			HaMat_double es,green;
	        mat_scale(es,nss,enr);
        	mat_diff(green,ng,es);
	        HaMat_double::mat_inverse(green);
	        q=19;
			for(m=1;m<=q;m++)
			{
			for(n=1;n<=q;n++)
			{
			for(a=1;a<=z;a++)
			{
			for(b=1;b<=z;b++)
			{
			tot1=tot1+(hh(nb-q*q+q*(n-1)+m,a)-enr*ss(nb-q*q+q*(n-1)+m,a))*green(a,b)*(hh(b,z+1)-enr*ss(b,z+1))*cos(2*pi*(k1*m+k2*n));
			tot2=tot2+(hh(nb-q*q+q*(n-1)+m,a)-enr*ss(nb-q*q+q*(n-1)+m,a))*green(a,b)*(hh(b,z+1)-enr*ss(b,z+1))*sin(2*pi*(k1*m+k2*n));
			}
			}
			}
			}
			tot=2*sin(2*pi*k3)*sin(2*pi*k3)*(tot1*tot1+tot2*tot2);
            tmat_el = tot;
	 cout << tot << endl;
	return True;
}
