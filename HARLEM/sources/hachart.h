/*! \file hachart.h

    2D and 3D Chart Classes in HARLEM based on PLPlot package
 
    \author Igor Kurnikov  
    \date 2011-
*/

#if !defined(HACHART_H)
#define HACHART_H


//#define WITH_PLPLOT 1 
#define WITH_PLPLOT 0 
#if(_DEBUG)
#define WITH_PLPLOT 0
#endif

#include "halinalg.h"

#if WITH_PLPLOT
class plstream;
#endif

class HaChart
//! Base Class for Chart classes 
{
public:
	HaChart();
	virtual ~HaChart();

	const char* GetTitle(); //!< Get Title of the Chart

#if WITH_PLPLOT
	virtual int ToPLpStream( plstream *pls ) = 0;  //!< Plot to PLplot stream
#endif

	virtual void Clear();     //!< Clear All Data and labeling info
	virtual void ClearData(); //!< Clear Data Arrays

protected:
	std::string title;  //!< Chart title

	double xmin_vp; //!< Viewport xmin coordinate in reduced subpage coordinate (from 0.0 to 1.0)
	double xmax_vp; //!< Viewport xmax coordinate in reduced subpage coordinate (from 0.0 to 1.0)
	double ymin_vp; //!< Viewport ymin coordinate in reduced subpage coordinate (from 0.0 to 1.0)
	double ymax_vp; //!< Viewport ymax coordinate in reduced subpage coordinate (from 0.0 to 1.0)
};

class HaChartAxis
//!< Base Chart Axis class
{
public:
	HaChartAxis();
	virtual ~HaChartAxis();

	bool pos_crd_orig;     //!< Draw axis at coordinate origin ( y = 0 for X-axis, x = 0 for Y-axis )
	bool pos_bot_left;     //!< Draw axis at frame bottom (for X-axis) or left edge of the frame ( Y or Z axis )
	bool pos_top_right;    //!< Draw axis at frame top (for X-axis) or right edge if the frame ( Y or Z axis )
	bool edge_line;        //!< Draw edge line if the axis is drawn at the frame edge
	bool ticks_invert;     //!< Inverts tick marks, so they are drawn outwards, rather than inwards
	bool ticks_major;      //!< Draw major ticks
	bool ticks_minor;      //!< Draw minor ticks
	bool grid_ticks_major; //!< Draw a grid at the major tick intervals
	bool grid_ticks_minor; //!< Draw a grid at the minor tick intervals
	bool lbl_major_ticks;  //!< Write labels at the major tick intervals 
	bool lbl_unconv_loc;   //!< Write labels at unconventional location ( above box edge for X, right of box for Y)
	bool lbl_custom;       //!< Use custom labeling function
	bool lbl_fixed_point;  //!< Always use fixed point numeric labels
	bool lbl_log;          //!< Labels axis logarithmically
	bool lbl_par_base;     //!< Write numeric labels parallel to the base of the graph rather than parallel to the axis (for Y axis) 

	std::string axis_label;  //!< label of the axis
	
	double ticks_major_interval;     //!< Interval between major ticks ( if =0.0 set automatically)
	int    ticks_minor_num_interval; //!< Number of subintervals between major X-axis ticks for minor ticks. ( if =0 set automatically)
};

class HaChart2D: public HaChart
//! 2D Chart class 
{
public:
	HaChart2D();
	virtual ~HaChart2D();

	virtual void Clear();     //!< Clear All Data and labeling info
	virtual void ClearData(); //!< Clear Data Arrays

#if WITH_PLPLOT
	virtual int ToPLpStream( plstream *pls );  //!< Plot to PLplot stream
#endif

	int AddXYData(const HaVec_double& x,const HaVec_double& y); //!< Add X-Y Data Series to the chart
	int AddYData(const HaVec_double& y); //!< Add an array of Y values (X is set to 0-based index)

	void SetAxesLbls( const std::string& xlbl, const std::string& ylbl); //!< Set Axes labels
	void SetAxisXLbl( const std::string& xlbl); //!< Set X axis label
	void SetAxisYLbl( const std::string& ylbl); //!< Set Y axis label

	std::string GetXLbl() const; //!< Get X axis label
	std::string GetYLbl() const; //!< Get Y axis label

	void GetXRange(double& xmin_n, double& xmax_n); //!< Get xmin and xmax for the data sets
	void GetYRange(double& ymin_n, double& ymax_n); //!< Get ymin and ymax for the data sets

	int GetNSer() const; //!< Get the number of X-Y data series

	HaVec_double* GetXData(int i_ser); //!< Get X data for series (i_ser) return NULL at error 
	HaVec_double* GetYData(int i_ser); //!< Get Y data for series (i_ser) return NULL at error 

protected:
	std::vector<HaVec_double> xc;
	std::vector<HaVec_double> yc;

	double xmin;
	double xmax;
	double ymin;
	double ymax;

	int xrange_auto_flag; //!< Flag to set xmin and xmax automatically for data sets 
	int yrange_auto_flag; //!< Flag to set ymin and ymax automatically for data sets

	HaChartAxis x_axis; //!< X-axis parameters
	HaChartAxis y_axis; //!< Y-axis parameters
};

class HaChart3D: public HaChart
//! 3D Chart class 
{
public:
	HaChart3D();
	virtual ~HaChart3D();

	void Clear();     //!< Clear All Data and labeling info
	void ClearData(); //!< Clear Data Arrays

#if WITH_PLPLOT
	virtual int ToPLpStream( plstream *pls );  //!< Plot to PLplot stream
#endif

	int AddXYZData(const HaVec_double& x,const HaVec_double& y, const HaVec_double& z); //!< Add X-Y-Z Data Series to the chart

	void SetAxesLbls( const std::string& xlbl, const std::string& ylbl,const std::string& zlbl); //!< Set Axes labels
	void SetAxisXLbl( const std::string& xlbl); //!< Set X axis label
	void SetAxisYLbl( const std::string& ylbl); //!< Set Y axis label
	void SetAxisZLbl( const std::string& ylbl); //!< Set Z axis label

	std::string GetXLbl() const; //!< Get X axis label
	std::string GetYLbl() const; //!< Get Y axis label
	std::string GetZLbl() const; //!< Get Z axis label

	void GetXRange(double& xmin_n, double& xmax_n); //!< Get xmin and xmax for the data sets
	void GetYRange(double& ymin_n, double& ymax_n); //!< Get ymin and ymax for the data sets
	void GetZRange(double& zmin_n, double& zmax_n); //!< Get ymin and ymax for the data sets

	int GetNSer() const; //!< Get the number of X-Y-Z data series

	HaVec_double* GetXData(int i_ser); //!< Get X data for series (i_ser) return NULL at error 
	HaVec_double* GetYData(int i_ser); //!< Get Y data for series (i_ser) return NULL at error 
	HaVec_double* GetZData(int i_ser); //!< Get Z data for series (i_ser) return NULL at error 

protected:
	std::vector<HaVec_double> xc;
	std::vector<HaVec_double> yc;
	std::vector<HaVec_double> zc;

	double xmin; //!< minimal X coordinate of the plotting box in world coordinates ( data coordinates) 
	double xmax; //!< maximal X coordinate of the plotting box in world coordinates ( data coordinates) 
	double ymin; //!< minimal Y coordinate of the plotting box in world coordinates ( data coordinates) 
	double ymax; //!< maximal Y coordinate of the plotting box in world coordinates ( data coordinates) 
	double zmin; //!< minimal Z coordinate of the plotting box in world coordinates ( data coordinates) 
	double zmax; //!< maximal Z coordinate of the plotting box in world coordinates ( data coordinates) 

	double xmin_win; //!< Viewing Window xmin coordinate (normalized coordinates)  
	double xmax_win; //!< Viewing Window xmax coordinate (normalized coordinates) 
	double ymin_win; //!< Viewing Window ymin coordinate (normalized coordinates) 
	double ymax_win; //!< Viewing Window ymax coordinate (normalized coordinates) 

	double xsize_box; //!< size of the plotting box along X axis (normalized coordinates)
	double ysize_box; //!< size of the plotting box along Y axis (normalized coordinates)
	double zsize_box; //!< size of the plotting box along Z axis (normalized coordinates)

	double alt_view_ang; //!<  altitude angle for location of the observer (viewing angle above the xy plane) (in degrees from 0.0 to 90.0 )
    double az_view_ang;  //!<  azimuth angle for location of the observer  (in degrees from 0.0 to 90.0 )
	                     //!< The azimuth is defined so that when az = 0, the observer sees the xz plane face on, 
	                     //!< and as the angle is increased, the observer moves clockwise around the box 
	                     //!< as viewed from above the xy plane. The azimuth can take on any value

	int xrange_auto_flag; //!< Flag to set xmin and xmax automatically for data sets 
	int yrange_auto_flag; //!< Flag to set ymin and ymax automatically for data sets
	int zrange_auto_flag; //!< Flag to set ymin and ymax automatically for data sets

	HaChartAxis x_axis; //!< X-axis parameters
	HaChartAxis y_axis; //!< Y-axis parameters
	HaChartAxis z_axis; //!< Z-axis parameters

};



class HaChartPanel
//! Class for a panel of 2D and 3D charts
{
public:
	HaChartPanel();
	virtual ~HaChartPanel();

#if WITH_PLPLOT
	int ToPLpStream( plstream *pls );  //!< Plot to PLplot stream
#endif
	
	HaChart2D* AddChart2D();    //!< Add new 2D chart to the panel
	HaChart3D* AddChart3D();    //!< Add new 3D chart to the panel
	HaChart* GetChart(int idx); //!< Get Chart by index  

protected:
	std::vector< HaChart* > chart_arr;
	
};

#endif  // !defined HACHART_H
