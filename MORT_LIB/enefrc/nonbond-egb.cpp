#include <common.hpp>
#include <object.hpp>
#include "nonbond.hpp"
#include "ctrlparm.hpp"

namespace mort
{
    using namespace boost;

    void get_reff(int natom, const double* x,const double* brnr,const double* gbfs,double* reff, const ctrlparm_t& p)
    {
        static const double ta = 1.0/3.0;
        static const double tb = 2.0/5.0;
        static const double tc = 3.0/7.0;
        static const double td = 4.0/9.0;
        static const double tdd= 5.0/11.0;
        

        for(int i=0; i < natom; ++i)
	{
	    reff[i] = 0.0;
	}

        double fsmax = 0.0;
	for(int i=0; i < natom; ++i)
        {
	    double ri = brnr[i]-p.offset;
            double si = ri*gbfs[i];
            if( fsmax < si ) fsmax = si;
	}


        double rgbmax1i = 1.0/p.rgbmax;
	double rgbmax2i = rgbmax1i * rgbmax1i;
	double rgbmaxpsmax2 = (p.rgbmax+fsmax)*(p.rgbmax+fsmax);

	for(int i=0; i < natom; ++i)
	{
	    double ri = brnr[i]-p.offset;
            double si = ri*gbfs[i];
	    double si2 = si*si;
            double xi = x[3*i];
	    double yi = x[3*i+1];
	    double zi = x[3*i+2];

            for(int j=i+1; j < natom; ++j)
	    {
                double xij = xi - x[3*j];
                double yij = yi - x[3*j+1];
		double zij = zi - x[3*j+2];
                double r2  = xij*xij + yij*yij + zij*zij;

                if( r2 > rgbmaxpsmax2 )
		    continue;

                double dij = sqrt(r2);
                double rj = brnr[j]-p.offset;
	        double sj = gbfs[j]*rj;
                double sj2 = sj*sj;
		
                double temp3 = 0.0;
                double temp4 = 0.0;

                if( dij <= p.rgbmax+sj )
		{
                    if( dij > p.rgbmax - sj )
                    {
                        double uij = 1.0/(dij-sj);
                        temp3 = temp3 - 0.125/dij * (1.0+2.0*dij*uij+rgbmax2i*(r2 -4.0*p.rgbmax*dij -sj2)+2.0*log((dij-sj)*rgbmax1i));
                    }
                    else if( dij > 4.0*sj )
                    {
                        double tmpsd = sj2/r2;
                        double dumbo = ta + tmpsd*(tb+tmpsd*(tc+tmpsd*(td+tmpsd*tdd)));
                        temp3 = temp3 - tmpsd*sj*dumbo/r2;                                
                    }
                    else
                    {
                        double tmp2 = dij + sj;
                        double tmp4 = 1.0;
                        
                        if( dij > ri+sj )
                        {
                            tmp4 = dij - sj;
                        }
                        else if( dij > std::abs(ri-sj) )
                        {
                            tmp4 = ri;
                        }
                        else if( ri < sj )
                        {
                            tmp4 = sj - dij;
                        }
                        else
                        {
                            tmp4 = 1.0;
                        }

                        tmp2 = 1.0/tmp2;
                        tmp4 = log( tmp2*tmp4 );

                        if( dij > ri+sj )
                        {
                            temp3 = temp3 - 0.5*( sj/(r2-sj2) + 0.5*tmp4/dij );    
                        }
                        else if( dij > std::abs(ri-sj) )
                        {
                            double theta = 0.5/ri/dij*(r2+ri*ri-sj2);
                            temp3 = temp3 - 0.25*( (2.0-theta)/ri - tmp2 + tmp4/dij );
                        }
                        else if( dij < sj )
                        {
                            temp3 = temp3 - 0.4*sj/(r2-sj2) + 1.0/ri + 0.25/dij*tmp4;
                        }
                        
                    }
                }
                
                if( dij <= p.rgbmax+si ) 
                {
                    if( dij > p.rgbmax - si )
                    {
                        double uij = 1.0/(dij-si);
                        temp4 = temp4 - 0.125/dij*(1.0 + 2.0 * dij*uij + rgbmax2i*(r2-4.0*p.rgbmax*dij-si2)+2.0*log( (dij-si)*rgbmax1i) );
                    }
                    else if( dij > 4.0*si )
                    {
                        double tmpsd = si2/r2;
                        double dumbo = ta + tmpsd*(tb+tmpsd*(tc+tmpsd*(td+tmpsd*tdd)));
                        temp4 = temp4 - tmpsd*si*dumbo/r2;
                    }
                    else
                    {
                        double tmp3 = dij + si;
                        double tmp5 = 1.0;
                        
                        if( dij > rj+si )
                        {
                            tmp5 = dij - si;
                        }
                        else if( dij > std::abs(rj-si) )
                        {
                            tmp5 = rj;
                        }
                        else if( rj < si ) 
                        {
                            tmp5 = si - dij;
                        }
                        else
                        {
                            tmp5 = 1.0;
                        }

                        tmp3 = 1.0/tmp3;
                        tmp5 = log( tmp3*tmp5 );
                        
                        if( dij>rj+si) 
                        {
                            temp4 = temp4 - 0.5*( si/(r2-si2) + 0.5/dij*tmp5 );
                        }
                        else if( dij > std::abs(rj-si) )
                        {
                            double theta = 0.5/rj/dij*(r2 + rj*rj - si2 );
                            temp4 = temp4 - 0.25*( (2.0-theta)/rj - tmp3 + tmp5/dij );
                        }
                        else if( rj < si ) 
                        {
                            temp4 = temp4 - 0.5*si/(r2-si2) + 1.0/rj + 0.25*tmp5/dij;
                        }
                    }
                    
                }

                reff[i] += temp3;
                reff[j] += temp4;
  	    }
        }

        for(int i=0; i<natom; ++i)
        {
            reff[i] += 1.0 / (brnr[i]-p.offset);
            reff[i] = std::max(1.0/reff[i], brnr[i]-p.offset);
        }
    }

        
    double get_egb(const nonbond_t& nb, const molecule_t& m, const ctrlparm_t& p)
    {
        assert( m.natom() == (int)nb.list.size() );

        vector<double> reff(m.natom());
	const vector<double>& vpos = get_vvec(m, ATOM, POSITION);
        const vector<double>& pchg = get_dvec(m, ATOM, PCHG);
        const vector<double>& brnr = get_dvec(m, ATOM, BRNR);
        const vector<double>& gbfs = get_dvec(m, ATOM, GBFS);
        get_reff(m.natom(), &vpos[0], &brnr[0], &gbfs[0], &reff[0], p);

	double egb = 0.0;
        int natom = m.natom();
        for(int i=0; i < natom; ++i)
        {
            double qi = pchg[i];
	    double bri = reff[i];
	    double qid = qi*DIELFAC;

	    for(int k=0; k<(int)nb.list[i].size(); ++k)
	    {
	        int j = nb.list[i][k];
                double qiqjd = qid*pchg[j];
                double bibj  = bri*reff[j];
                double rij2  = nb.dis2[i][k];
                double temp  = qiqjd / sqrt(rij2 + bibj * exp(-0.25*rij2/bibj) );
                egb = egb - temp;
            }
            double qid2h = 0.5*qid*qi;
            egb = egb - qid2h/bri;
        }
        return egb*INVCHG2;
    }

    double get_egb(const nonbond_t& nb, const atomvec_t& vec1, const atomvec_t& vec2, const ctrlparm_t& p)
    {
        assert( vec1.size()>0 && vec2.size()>0 );
            
        int natom = nb.list.size();
        vector<double> reff(natom);
	const vector<double>& vpos = get_vvec(vec1[0].getmol(), ATOM, POSITION);
        const vector<double>& pchg = get_dvec(vec1[0].getmol(), ATOM, PCHG);
        const vector<double>& brnr = get_dvec(vec1[0].getmol(), ATOM, BRNR);
	const vector<double>& gbfs = get_dvec(vec1[0].getmol(), ATOM, GBFS);
	get_reff(vec1[0].getmol().natom(), &vpos[0], &brnr[0], &gbfs[0], &reff[0], p);

        double egb=0.0;
	atomvec_t::const_iterator ai = vec1.begin();
	for( ; ai != vec1.end(); ++ai )
	{
	    int i = ai->absid();
	    atomvec_t::const_iterator aj = vec2.begin();
	    for( ; aj != vec2.end(); ++aj )
	    {
	        int j = aj->absid();

                int m = std::min(i, j);
                int n = std::max(i, j);

                vector<int>::const_iterator p = std::find(nb.list[m].begin(), nb.list[m].end(), n);
                if( p == nb.list[m].end() ) continue; 

	        int k = p - nb.list[m].begin();
                double qiqjd = pchg[m] * pchg[n] * DIELFAC;
                double bibj  = reff[m] * reff[n];
                double rij2  = nb.dis2[m][k];
                double temp  = qiqjd / sqrt(rij2 + bibj * exp(-0.25*rij2/bibj) );
                egb = egb - temp;
            }
        }

        return egb*INVCHG2;
    }

    numvec get_enb(const nonbond_t& nb, const molecule_t& m) 
    {
        const vector<double>& pchg = get_dvec(m, ATOM, PCHG);
        const vector<double>& vdwr = get_dvec(m, ATOM, RSTAR);
        const vector<double>& well = get_dvec(m, ATOM, DEPTH);

        double eelt=0.0;
        double evdw=0.0;
        int natom=nb.list.size();
        for(int i=0; i < natom;++i)
        {
            for(int k=0; k < (int)nb.list[i].size(); ++k)
      	    {
	        if( nb.type[i][k]==0 )
	        {
	            continue;
	        }

	        int j = nb.list[i][k];
	        double r1 = vdwr[i];
	        double r2 = vdwr[j];
                double w1 = well[i];
	        double w2 = well[j];
	        double rsq = (r1+r2)*(r1+r2);
	        double wsq = w1 * w2;
	        double f2 = rsq/nb.dis2[i][k];
	        double f6 = f2*f2*f2;
	        eelt += INVCHG2*pchg[i]*pchg[j]/nb.dist[i][k];
                evdw += sqrt(wsq)*f6*(f6-2.0);
	    }
	}

	return makevec(eelt,evdw);
    }

    numvec get_enb(const nonbond_t& nb, const atomvec_t& vec1, const atomvec_t& vec2)
    {
        assert( vec1.size()>0 && vec2.size()>0 );
        const vector<double>& pchg = get_dvec(vec1[0].getmol(), ATOM, PCHG);
        const vector<double>& vdwr = get_dvec(vec1[0].getmol(), ATOM, RSTAR);
        const vector<double>& well = get_dvec(vec1[0].getmol(), ATOM, DEPTH);
        double eelt = 0.0;
        double evdw = 0.0;
        atomvec_t::const_iterator ai = vec1.begin();
	for( ; ai != vec1.end(); ++ai )
	{
	    int i = ai->absid();
	    atomvec_t::const_iterator aj = vec2.begin();
	    for( ; aj != vec2.end(); ++aj)
	    {
	        int j = aj->absid();

                int m = std::min( i, j );
                int n = std::max( i, j );

                vector<int>::const_iterator p=std::find(nb.list[m].begin(), nb.list[m].end(), n);
	        if( p == nb.list[m].end() ) continue;

                int k = p - nb.list[m].begin();
	        if( nb.type[m][k] == 0 ) continue;

	        double r1 = vdwr[m];
	        double r2 = vdwr[n];
                double w1 = well[m];
                double w2 = well[n];
	        double rsq = (r1+r2)*(r1+r2);
	        double wsq = w1 * w2;
	        double f2 = rsq/nb.dis2[m][k];
	        double f6 = f2*f2*f2;
	        eelt += INVCHG2*pchg[m]*pchg[n]/nb.dist[m][k];
                evdw += sqrt(wsq)*f6*(f6-2.0);
            }
        }

        return makevec(eelt, evdw);
    }          
    
    numvec nonbond_egb(const molecule_t& m, const ctrlparm_t& p)
    {
        nonbond_t nb(m);

        nb.list_for_egb(p.cut);

        double egb = get_egb(nb, m, p);

        numvec enb = get_enb(nb, m);

        return makevec(enb[0], enb[1], egb);
    }

    numvec nonbond_egb(const atomvec_t& vec1, const atomvec_t& vec2, const ctrlparm_t& p)
    {
        assert( vec1.size()>0 && vec2.size()>0 );

        nonbond_t nb( vec1[0].getmol() );
	nb.list_for_egb(p.cut);

	double egb = get_egb(nb, vec1, vec2, p);
        numvec enb = get_enb(nb, vec1, vec2);

        return makevec(enb[0], enb[1], egb);
    }

} // namespace mort
