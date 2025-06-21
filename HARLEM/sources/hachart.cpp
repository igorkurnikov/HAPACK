/*! \file hachart.cpp

    2D and 3D Plotting Classes in HARLEM Based on PLPlot package
 
    \author Igor Kurnikov  
    \date 2011-
*/

#include "haio.h"
#include <stdexcept>
#include <memory>


// common includes

#include "halinalg.h"
#include "hachart.h"

#if WITH_PLPLOT
#include "plplot/plstream.h"
#include "plplot/plevent.h"
#include "plplot/wxPLplotstream.h"
#endif

HaChart::HaChart()
{
	Clear();
}

HaChart::~HaChart()
{

}

const char* HaChart::GetTitle()
{
	return title.c_str();
}

void HaChart::Clear()
{
	xmin_vp = 0.0;
	xmax_vp = 0.0; 
	ymin_vp = 0.0; 
	ymax_vp = 0.0;

	title.clear();
}

void HaChart::ClearData()
{

}

HaChartAxis::HaChartAxis()
{
	pos_crd_orig = false;  
	pos_bot_left   = true;    
	pos_top_right  = true;  
    edge_line = true;        
    ticks_invert = false;     
    ticks_major  = true;      
    ticks_minor  = true;     
    grid_ticks_major = false; 
    grid_ticks_minor = false;
    lbl_major_ticks = true; 
    lbl_unconv_loc  = false;   
    lbl_custom = false;   
	lbl_fixed_point = false;
    lbl_log = false;  
	lbl_par_base = true;

	axis_label = "Data"; 
	
    ticks_major_interval = 0.0;
    ticks_minor_num_interval  = 0;
}
	
HaChartAxis::~HaChartAxis()
{

}


HaChart2D::HaChart2D()
{
	Clear();
}

HaChart2D::~HaChart2D()
{
	
}

void HaChart2D::Clear()
{
	HaChart::Clear();
	
	x_axis.axis_label = "x";
	y_axis.axis_label = "y";

	xrange_auto_flag = TRUE;
	yrange_auto_flag = TRUE;

	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;

	ClearData();
}
	
void HaChart2D::ClearData()
{
	xc.clear();
	yc.clear();
}

#if WITH_PLPLOT
int HaChart2D::ToPLpStream( plstream *pls )
{
	double xmin,xmax,ymin,ymax;
	GetXRange(xmin,xmax);
	GetYRange(ymin,ymax);
			
	std::string title = GetTitle();
	std::string x_lbl = GetXLbl();
	std::string y_lbl = GetYLbl();

	std::string xopt;
	std::string yopt;

// Set X axis parameters

	if(x_axis.pos_crd_orig)    xopt.append("a");
	if(x_axis.pos_bot_left)    
	{
		if( x_axis.edge_line ) xopt.append("b");
		else xopt.append("u");
	}
	if(x_axis.pos_top_right)   
	{
		if( x_axis.edge_line ) xopt.append("c");
		else xopt.append("w");
	}
	if(x_axis.lbl_fixed_point) xopt.append("f");
	if(x_axis.grid_ticks_major) xopt.append("g");
	if(x_axis.grid_ticks_minor) xopt.append("h");
	if(x_axis.ticks_invert)     xopt.append("i");
	if(x_axis.lbl_log)          xopt.append("l");
	if(x_axis.lbl_major_ticks)
	{
		if( x_axis.lbl_unconv_loc ) xopt.append("m");
		else xopt.append("n");
	}
//	if(x_axis.lbl_custom) xopt.append("o"); // Custom labeling are not set yet...
	if(x_axis.ticks_major) xopt.append("t");
	if(x_axis.ticks_minor) xopt.append("s");

// Set Y axis parameters

	if(y_axis.pos_crd_orig)    yopt.append("a");
	if(y_axis.pos_bot_left)    
	{
		if( y_axis.edge_line ) yopt.append("b");
		else yopt.append("u");
	}
	if(y_axis.pos_top_right)   
	{
		if( y_axis.edge_line ) yopt.append("c");
		else yopt.append("w");
	}
	if(y_axis.lbl_fixed_point)  yopt.append("f");
	if(y_axis.grid_ticks_major) yopt.append("g");
	if(y_axis.grid_ticks_minor) yopt.append("h");
	if(y_axis.ticks_invert)     yopt.append("i");
	if(y_axis.lbl_log)          yopt.append("l");
	if(y_axis.lbl_major_ticks)
	{
		if( y_axis.lbl_unconv_loc ) yopt.append("m");
		else yopt.append("n");
	}
//	if(y_axis.lbl_custom) yopt.append("o"); // Custom labeling are not set yet...
	if(y_axis.ticks_major)   yopt.append("t");
	if(y_axis.ticks_minor)   yopt.append("s");
	if(y_axis.lbl_par_base ) yopt.append("v");


	double xtick = x_axis.ticks_major_interval;
	double ytick = y_axis.ticks_major_interval;
	int    nxsub = x_axis.ticks_minor_num_interval;
	int    nysub = y_axis.ticks_minor_num_interval;

	pls->col0( 1 );
//	pls->env( xmin, xmax, ymin, ymax, 0, 0 );
	pls->adv(0);
	if( ( xmin_vp < (xmax_vp -0.001) ) && ( ymin_vp < (ymax_vp -0.001) ) )
	{
		pls->vpor(xmin_vp,xmax_vp,ymin_vp,ymax_vp);
	}
	else
	{
		pls->vsta(); // Standard Viewport
	}
	pls->wind(xmin,xmax,ymin,ymax);
	pls->box(xopt.c_str(),xtick,nxsub,yopt.c_str(),ytick,nysub); 

	pls->col0( 2 );
	pls->lab( x_lbl.c_str(), y_lbl.c_str(), title.c_str() );

	int ns = GetNSer();
	int is;
	for(is = 0; is < ns; is++)
	{
		int icol = is % 12 + 2;
		pls->col0( icol );

		HaVec_double* px = GetXData(is);
		HaVec_double* py = GetYData(is);

		int np = py->size();

		if( px->size() != np )
		{
			PrintLog("X and Y arrays has different dimensions \n");
			continue;
		}

		pls->width( 2 );
		pls->line( np, px->v(), py->v());
	}
	return TRUE;
}
#endif
	
int HaChart2D::AddXYData(const HaVec_double& x,const HaVec_double& y)
{
	try
	{
		if(y.size() == 0) throw std::runtime_error(" array y is empty ");
		if(x.size() != y.size() ) throw std::runtime_error(" Dimensions of x and y arrays are not equal ");
		
		xc.push_back(x);
		yc.push_back(y);
	}
	catch( std::exception& ex)
	{
		PrintLog("Error in HaChart2D::AddXYdata() \n");
		PrintLog("%s\n",ex.what());

		return FALSE;
	}
	return TRUE;
}

int HaChart2D::AddYData(const HaVec_double& y)
{
	try
	{
		if(y.size() == 0) throw std::runtime_error(" array y is empty ");
		int np = y.size();
		HaVec_double x(np);
		int i;
		for(i = 0; i < np; i++)
		{
			x[i] = (double) i;
		}
		xc.push_back(x);
		yc.push_back(y);
	}
	catch( std::exception& ex)
	{
		PrintLog("Error in HaChart2D::AddYdata() \n");
		PrintLog("%s\n",ex.what());

		return FALSE;
	}
	return TRUE;
}	

void HaChart2D::SetAxesLbls( const std::string& xlbl, const std::string& ylbl)
{
	SetAxisXLbl(xlbl);
	SetAxisYLbl(ylbl);
}

void HaChart2D::SetAxisXLbl( const std::string& xlbl)
{
	x_axis.axis_label = xlbl;
}
	

void HaChart2D::SetAxisYLbl( const std::string& ylbl)
{ 
	y_axis.axis_label = ylbl;
}


void HaChart2D::GetXRange(double& xmin_n, double& xmax_n)
{
	const double x_small = -1.0e12;
	const double x_large =  1.0e12;

	if( xrange_auto_flag || (xmax - xmin) < 1.0e-12 )
	{
		xmin = x_large;
		xmax = x_small;
		int ns = xc.size();
		int is;
		for( is = 0; is < ns; is++)
		{
			int n = xc[is].size();
			int i;
			for( i = 0; i < n; i++ )
			{
				if( xc[is][i] <  xmin ) xmin = xc[is][i];
				if( xc[is][i] >  xmax ) xmax = xc[is][i];
			}
		}
		if( xmax < xmin )
		{
			xmin = 0.0;
			xmax = 1.0;
		}
	}
	
	xmin_n = xmin;
	xmax_n = xmax;
}

void HaChart2D::GetYRange(double& ymin_n, double& ymax_n)
{
	const double y_small = -1.0e12;
	const double y_large =  1.0e12;

	if( yrange_auto_flag || (ymax - ymin) < 1.0e-12 )
	{
		ymin = y_large;
		ymax = y_small;
		int ns = yc.size();
		int is;
		for( is = 0; is < ns; is++)
		{
			int n = yc[is].size();
			int i;
			for( i = 0; i < n; i++ )
			{
				if( yc[is][i] <  ymin ) ymin = yc[is][i];
				if( yc[is][i] >  ymax ) ymax = yc[is][i];
			}
		}
		if( ymax < ymin )
		{
			ymin = 0.0;
			ymax = 1.0;
		}
	}
	
	ymin_n = ymin;
	ymax_n = ymax;
}

std::string HaChart2D::GetXLbl() const
{
	return x_axis.axis_label;
}

std::string HaChart2D::GetYLbl() const
{
	return  y_axis.axis_label;
}

int HaChart2D::GetNSer() const
{
	return yc.size();
}

HaVec_double* HaChart2D::GetXData(int i_ser)
{
	int ns = GetNSer();
	if( i_ser < 0 || i_ser >= ns ) return NULL;
	if( xc.size() == 0 ) return NULL;
	if( i_ser >= xc.size() ) return &xc[0];
	return &xc[i_ser];
}

HaVec_double* HaChart2D::GetYData(int i_ser)
{
	int ns = GetNSer();
	if( i_ser < 0 || i_ser >= ns ) return NULL;
	return &yc[i_ser];
}


HaChart3D::HaChart3D()
{
	Clear();
}

HaChart3D::~HaChart3D()
{
	
}

void HaChart3D::Clear()
{
	HaChart::Clear();
	
	x_axis.axis_label = "x";
	y_axis.axis_label = "y";
	z_axis.axis_label = "z";

	xrange_auto_flag = TRUE;
	yrange_auto_flag = TRUE;
	zrange_auto_flag = TRUE;

	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;
	zmin = 0.0;
	zmax = 0.0;

	xmin_win = -2.5; 
	xmax_win =  2.5; 
	ymin_win = -2.5; 
	ymax_win =  4.0;

	xsize_box = 2.0; 
	ysize_box = 4.0; 
    zsize_box = 3.0; 

	alt_view_ang = 45.0; 
    az_view_ang  = 30.0;

	ClearData();
}
	
void HaChart3D::ClearData()
{
	xc.clear();
	yc.clear();
	zc.clear();
}

#if WITH_PLPLOT
int HaChart3D::ToPLpStream( plstream *pls )
{
	double xmin,xmax,ymin,ymax,zmin,zmax;
	GetXRange(xmin,xmax);
	GetYRange(ymin,ymax);
	GetZRange(zmin,zmax);
			
	std::string title = GetTitle();
	std::string x_lbl = GetXLbl();
	std::string y_lbl = GetYLbl();
	std::string z_lbl = GetZLbl();

	double alt = 45.0;
	double az  = 30.0;

	pls->adv(0);
	pls->vpor(xmin_vp,xmax_vp,ymin_vp,ymax_vp);
	pls->wind(xmin_win,xmax_win,ymin_win,ymax_win);
	pls->col0( 1 );
	
	pls->w3d(xsize_box,ysize_box,zsize_box,
		     xmin,xmax,ymin,ymax,zmin,zmax,
			 alt_view_ang,az_view_ang);  

	pls->box3("bnstu", x_lbl.c_str(), 0.0, 0,
		      "bnstu", y_lbl.c_str(), 0.0, 0,
              "bcdmnstuv", z_lbl.c_str(), 0.0, 0);

	int ns =  GetNSer();
	int is;
	for(is = 0; is < ns; is++)
	{
		int icol = is % 12 + 3;
		pls->col0( icol );

		HaVec_double* px = GetXData(is);
		HaVec_double* py = GetYData(is);
		HaVec_double* pz = GetZData(is);

		int np = pz->size();

		if( px->size() != np || py->size() != np )
		{
			PrintLog("X, Y, Z arrays have different dimensions \n");
			continue;
		}

		pls->width( 2 );
		pls->line3( np, px->v(), py->v(), pz->v() );
	}
	return TRUE;
}
#endif
	
int HaChart3D::AddXYZData(const HaVec_double& x,const HaVec_double& y, const HaVec_double& z)
{
	try
	{
		if(y.size() == 0) throw std::runtime_error(" array y is empty ");
		if(z.size() == 0) throw std::runtime_error(" array z is empty ");
		if(x.size() != y.size() || x.size() != z.size()) throw std::runtime_error(" Dimensions of x,y and z arrays are not equal ");
		
		xc.push_back(x);
		yc.push_back(y);
		zc.push_back(z);
	}
	catch( std::exception& ex)
	{
		PrintLog("Error in HaChart3D::AddXYZdata() \n");
		PrintLog("%s\n",ex.what());

		return FALSE;
	}
	return TRUE;
}

void HaChart3D::SetAxesLbls( const std::string& xlbl, const std::string& ylbl, const std::string& zlbl)
{
	SetAxisXLbl(xlbl);
	SetAxisYLbl(ylbl);
	SetAxisZLbl(zlbl);
}

void HaChart3D::SetAxisXLbl( const std::string& xlbl)
{
	x_axis.axis_label = xlbl;
}
	

void HaChart3D::SetAxisYLbl( const std::string& ylbl)
{
	y_axis.axis_label = ylbl;
}

void HaChart3D::SetAxisZLbl( const std::string& zlbl)
{
	z_axis.axis_label = zlbl;
}

void HaChart3D::GetXRange(double& xmin_n, double& xmax_n)
{
	const double x_small = -1.0e12;
	const double x_large =  1.0e12;

	if( xrange_auto_flag || (xmax - xmin) < 1.0e-12 )
	{
		xmin = x_large;
		xmax = x_small;
		int ns = xc.size();
		int is;
		for( is = 0; is < ns; is++)
		{
			int n = xc[is].size();
			int i;
			for( i = 0; i < n; i++ )
			{
				if( xc[is][i] <  xmin ) xmin = xc[is][i];
				if( xc[is][i] >  xmax ) xmax = xc[is][i];
			}
		}
		if( xmax < xmin )
		{
			xmin = 0.0;
			xmax = 1.0;
		}
	}
	
	xmin_n = xmin;
	xmax_n = xmax;
}

void HaChart3D::GetYRange(double& ymin_n, double& ymax_n)
{
	const double y_small = -1.0e12;
	const double y_large =  1.0e12;

	if( yrange_auto_flag || (ymax - ymin) < 1.0e-12 )
	{
		ymin = y_large;
		ymax = y_small;
		int ns = yc.size();
		int is;
		for( is = 0; is < ns; is++)
		{
			int n = yc[is].size();
			int i;
			for( i = 0; i < n; i++ )
			{
				if( yc[is][i] <  ymin ) ymin = yc[is][i];
				if( yc[is][i] >  ymax ) ymax = yc[is][i];
			}
		}
		if( ymax < ymin )
		{
			ymin = 0.0;
			ymax = 1.0;
		}
	}
	
	ymin_n = ymin;
	ymax_n = ymax;
}

void HaChart3D::GetZRange(double& zmin_n, double& zmax_n)
{
	const double z_small = -1.0e12;
	const double z_large =  1.0e12;

	if( zrange_auto_flag || (zmax - zmin) < 1.0e-12 )
	{
		zmin = z_large;
		zmax = z_small;
		int ns = zc.size();
		int is;
		for( is = 0; is < ns; is++)
		{
			int n = zc[is].size();
			int i;
			for( i = 0; i < n; i++ )
			{
				if( zc[is][i] <  zmin ) zmin = zc[is][i];
				if( zc[is][i] >  zmax ) zmax = zc[is][i];
			}
		}
		if( zmax < zmin )
		{
			zmin = 0.0;
			zmax = 1.0;
		}
	}
	
	zmin_n = zmin;
	zmax_n = zmax;
}


std::string HaChart3D::GetXLbl() const
{
	return x_axis.axis_label;
}

std::string HaChart3D::GetYLbl() const
{
	return y_axis.axis_label;
}

std::string HaChart3D::GetZLbl() const
{
	return z_axis.axis_label;
}

int HaChart3D::GetNSer() const
{
	return zc.size();
}

HaVec_double* HaChart3D::GetXData(int i_ser)
{
	int ns = GetNSer();
	if( i_ser < 0 || i_ser >= ns ) return NULL;
	return &xc[i_ser];
}

HaVec_double* HaChart3D::GetYData(int i_ser)
{
	int ns = GetNSer();
	if( i_ser < 0 || i_ser >= ns ) return NULL;
	return &yc[i_ser];
}

HaVec_double* HaChart3D::GetZData(int i_ser)
{
	int ns = GetNSer();
	if( i_ser < 0 || i_ser >= ns ) return NULL;
	return &zc[i_ser];
}

HaChartPanel::HaChartPanel()
{

}

HaChartPanel::~HaChartPanel()
{
	int ns = chart_arr.size();
	int is;
	for(is = 0; is < ns; is++)
	{
		delete chart_arr[is];
	}
}

HaChart2D* HaChartPanel::AddChart2D()
{
	chart_arr.push_back( new HaChart2D() );
	return (HaChart2D*) chart_arr.back();
}

HaChart3D* HaChartPanel::AddChart3D()
{
	chart_arr.push_back( new HaChart3D() );
	return (HaChart3D*) chart_arr.back();
}


HaChart* HaChartPanel::GetChart(int idx)
{
	if( idx < 0 || idx >= chart_arr.size() ) return NULL;
	return chart_arr[idx];
}

#if WITH_PLPLOT
int HaChartPanel::ToPLpStream( plstream *pls )
{
	int np = chart_arr.size();
	int ip;
	for(ip = 0; ip < np; ip++)
	{
		chart_arr[ip]->ToPLpStream(pls);
	}
	return TRUE;
}
#endif 

