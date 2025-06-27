/*! \file colors.cpp

    Classes to manage colors

   \author Igor Kurnikov
   \date 2004-

*/

#include <boost/algorithm/string.hpp>

#include "canvas3d.h"
#include "hamolview.h"
#include "math_num.h"

int RegisterStdColorNames()
{
	HaColor::RegisterColorName("WHITE",255,255,255);
	HaColor::RegisterColorName("BLACK",0,0,0);
	HaColor::RegisterColorName("RED",255,0,0);
	HaColor::RegisterColorName("GREEN",0,128,0);
	HaColor::RegisterColorName("BLUE",0,0,255);
	HaColor::RegisterColorName("YELLOW",255,255,0);
	HaColor::RegisterColorName("MAGENTA",255,0,255);
	HaColor::RegisterColorName("CYAN",0,255,255);

	HaColor::RegisterColorName("ORANGE",255,165,0);
	HaColor::RegisterColorName("GOLD",255,215,0);
	HaColor::RegisterColorName("LIME",0,255,0);
	HaColor::RegisterColorName("CHOCOLATE",210,105,30);
	HaColor::RegisterColorName("PURPLE",128,0,128);
	HaColor::RegisterColorName("GRAY",128,128,128);
	HaColor::RegisterColorName("DARKGRAY",169,169,169);
	HaColor::RegisterColorName("LIGHTGRAY",211,211,211);
	HaColor::RegisterColorName("BROWN",165,42,42);
	HaColor::RegisterColorName("PINK",255,192,203);
	HaColor::RegisterColorName("CORAL",255,127,80);

	return TRUE;
}

std::vector<ColorVal> HaColor::used_colors;
IntIntMap HaColor::cval_idx_map;  
std::map<ColorVal,std::string>    HaColor::col_name_map; 
std::map<std::string, ColorVal> HaColor::name_col_map;

int HaColor::color_names_init_flag = RegisterStdColorNames();

HaColor HaMolView::BackColor(0,0,0);
HaColor HaMolView::LabelColor(255,255,255);
HaColor HaMolView::BoxColor(255,255,255);   

HaColor::HaColor(int r, int g, int b)
{
	SetColor(r,g,b);
}

HaColor::HaColor(ColorVal cval)
{
	int rn = RComp(cval);
	int gn = GComp(cval);
	int bn = BComp(cval);
	SetColor(rn,gn,bn);
}

HaColor::HaColor(const HaColor& ref)
{
	r = ref.r;
	g = ref.g;
	b = ref.b;
	cval = ref.cval;
	cidx = ref.cidx;
}

HaColor::~HaColor()
{

}

int HaColor::SetColor(int rn, int gn, int bn)
{
	r = rn; g = gn; b = bn;
	cidx = RegisterColor(r, g, b);
    cval = GetColorVal(r,g,b);
	return cidx;
}

int HaColor::SetColor( const std::string& col_str )
{
	if( !ColorNameExist(col_str ) ) return FALSE;
	ColorVal cval = GetColorVal( col_str );
	int rn = RComp(cval);
	int gn = GComp(cval);
	int bn = BComp(cval);
	SetColor(rn,gn,bn);
	return TRUE;
}

bool HaColor::ColorNameExist(const std::string& col_str_par )
{
	std::string col_str = boost::trim_copy(col_str_par);
	boost::to_upper(col_str);
	if( name_col_map.count(col_str) == 0 ) return false;
	return true;
}

ColorVal HaColor::GetColorVal( const std::string& col_str_par )
{
	ColorVal cval = 0;

	std::string col_str = boost::trim_copy(col_str_par);
	boost::to_upper(col_str);

	if( name_col_map.count(col_str) == 0 ) return cval;
	return name_col_map[col_str];
}


std::string HaColor::GetColorName( int r, int g, int b )
{
	ColorVal cval =  HaColor::GetColorVal(r,g,b);
	if( col_name_map.count(cval) == 0 ) return "";
	return col_name_map[cval];
}

bool HaColor::operator==(const HaColor& ref) const 
{ 
	return cval == ref.cval; 
} 

static double Power( double x, int y )
{
    double result;

    result = x;
    while( y>1 )
    {   if( y&1 ) { result *= x; y--; }
        else { result *= result; y>>=1; }
    }
    return( result );
}

int HaColor::RegisterColor(int r, int g, int b)
{
	ColorVal cval = GetColorVal(r,g,b);
	if(cval_idx_map.count(cval) != 0) 
	{
		return cval_idx_map[cval];
	}

    double diffuse,fade;
    double temp,inten;
    int col;
    int i,k;

//	PrintLog(" Register Color %d %d %d \n", r,g,b);

	diffuse = 1.0 - HaMolView::Ambient;
	col = Canvas3D::empty_lut_idx;

//	PrintLog(" Current Canvas3D::empty_lut_idx = %d \n", Canvas3D::empty_lut_idx);
	for( i=0; i < ColourDepth; i++ )
	{   
		temp = (double)i/ColourMask;
		inten = diffuse*temp + HaMolView::Ambient;  //!< Intensity of the color interpolating from Ambient to 1.0
		fade = 1.0-inten;
		
		if( HaMolView::FakeSpecular )
		{	   
			temp = Power(temp,HaMolView::SpecPower);
			k = (int)(255*temp);
			temp = 1.0 - temp;
			inten *= temp;
			fade *= temp;
		}

		int ri,gi,bi;
		
		if( HaMolView::UseBackFade )
		{   
			temp = 1.0-inten;
			ri = (int)(r*inten + fade*HaMolView::BackColor.r); 
			gi = (int)(g*inten + fade*HaMolView::BackColor.g);
			bi = (int)(b*inten + fade*HaMolView::BackColor.b);
		} 
		else
		{   
			ri = (int)(r*inten); 
			gi = (int)(g*inten);
			bi = (int)(b*inten);
		}
		
		if( HaMolView::FakeSpecular )
		{   
			ri += k;
			gi += k;
			bi += k;
		}
		
//		if(i == 0)
//		{
//			ri = r;
//			gi = g;
//			bi = b;
//		}

		ColorVal col_val = GetColorVal( ri, gi, bi );

//		PrintLog("Setting color (%d,%d,%d)= %d at index= %d \n",ri, gi, bi, col_val, col+i);

		Canvas3D::Lut[col+i] = GetColorVal( ri, gi, bi );
	}

	used_colors.push_back(cval);
	cval_idx_map[cval] = col; 

	Canvas3D::empty_lut_idx += ColourDepth;
	return col;
}

ColorVal HaColor::GetColorVal(int r, int g, int b)
{
	return ( (ColorVal)((r<<8)|g)<<8 ) | b;	
}

ColorVal HaColor::GetPureColorByCIdx(int cidx)
{
	int idx_pure_col = cidx/ColourDepth;
	return used_colors[idx_pure_col];
}

void HaColor::GetPureRGBbyCIdx(int cidx, int& r, int& g, int& b )
{
	int idx_pure_col = cidx/ColourDepth;
	ColorVal cval = used_colors[idx_pure_col];
	r = RComp(cval);
	g = GComp(cval);
	b = BComp(cval);
}

void HaColor::GetRGBFromColVal(ColorVal cval, int& r, int& g, int& b )
{
	r = RComp(cval);
	g = GComp(cval);
	b = BComp(cval);
}

int HaColor::RegisterColorName(const std::string& col_name_par, int r, int g, int b)
{
	std::string col_name = boost::trim_copy(col_name_par);
	boost::to_upper(col_name);
	ColorVal cval = GetColorVal(r,g,b);
	if( name_col_map.count(col_name) > 0 ) return FALSE;

	col_name_map[cval]     = col_name;
	name_col_map[col_name] = cval;
	
	return TRUE;
}

HaColorMap::HaColorMap()
{

}


HaColorMap::~HaColorMap()
{

}

int HaColorMap::AddColor(int r, int g, int b)
{
	HaColor color(r,g,b);
	colors.push_back(color);
	return (colors.size() - 1);
}


int HaColorMap::AddUniformRange(int num_cols, HaColor& col1, HaColor& col2)
{
	if(num_cols < 2) return FALSE;
	
	int i;
	for( i = 0; i < num_cols; i++)
	{
		if( i== 0)
		{
			if( colors.size() > 0 && (colors.back()).cval == col1.cval) 
				continue;  // if col1 is the same as the lat color in the sequence do not add
		}
		int r = (col1.r*(num_cols-i) + col2.r*i)/num_cols;
		int g = (col1.g*(num_cols-i) + col2.g*i)/num_cols;
		int b = (col1.b*(num_cols-i) + col2.b*i)/num_cols;
		HaColor color(r,g,b);
        AddColor(r,g,b);
	}
	return TRUE;
}


int StrColorMap::AddStrColorPair( const char* str, int r, int g, int b)
{
	std::string strc(str);
    ColorVal cval = HaColor::GetColorVal(r,g,b);
	int idx;
	if( cval_idx_map.count(cval) == 0)
	{
		idx = AddColor(r,g,b);	
        cval_idx_map[cval] = idx;
	}
	else
	{
		idx = cval_idx_map[cval];
	}
	str_idx_map[strc] = idx;
	return idx;
}

HaColor* StrColorMap::GetColorForStr(const char* str)
{
	if( str_idx_map.count(str) == 0)
	{
		return NULL;
	}
	int idx = str_idx_map[str];

	return &colors[idx];
}

IValColorMap::IValColorMap()
{

}

IValColorMap::~IValColorMap()
{

}


int IValColorMap::AddIntColorPair( int ival, int r, int g, int b)
{
    ColorVal cval = HaColor::GetColorVal(r,g,b);
	int idx;
	if( cval_idx_map.count(ival) == 0)
	{
		idx = AddColor(r,g,b);	
        cval_idx_map[cval] = idx;
	}
	else
	{
		idx = cval_idx_map[cval];
	}
	int_idx_map[ival] = idx;
	return idx;
}

HaColor* IValColorMap::GetColorForInt(int ival)
{
	if( int_idx_map.count(ival) == 0)
	{
		return NULL;
	}
	int idx = int_idx_map[ival];

	return &colors[idx];
}


DRangeColorMap::DRangeColorMap()
{
	min_val = 0.0;
	max_val = 0.0;
}

DRangeColorMap::~DRangeColorMap()
{

}

HaColor* DRangeColorMap::GetColorForVal(double val)
{
	int n = GetNColors();
	if( n == 0) return NULL;

	if( val < min_val) return &colors[0];
	if( val > max_val) return &colors[n-1];
	if( max_val < (min_val + DBL_EPSILON)) return &colors[0];

	int idx = nintf( (n-1)*(val - min_val)/(max_val - min_val) );
	return &colors[idx];
}


DValColorMap::DValColorMap()
{
	Clear();
}

DValColorMap::~DValColorMap()
{

}

void DValColorMap::Clear()
{
	colors.clear();
	min_values.clear();
}

bool DValColorMap::IsValid() const
{
	if( colors.size() == 0 ) return false;
	if( min_values.size() != colors.size() ) return false;
	return true;
}

int DValColorMap::AddColorAndMinVal(int r, int g, int b, double val_min, int n_interp_col )
{
	int ires;
	try
	{
		if( colors.size() == 0 ) val_min = -1.0e10;

		int nl = min_values.size();
		
		if( nl > 0 && min_values[nl-1] > val_min ) throw std::runtime_error("Low limit of the added color is smaller than limit of the previous color");
		
		if( n_interp_col > 0  && colors.size() >= nl && nl > 0)
		{
			double prev_val_min = min_values[nl-1];
			double prev_r = colors[nl-1].r;
			double prev_g = colors[nl-1].g;
			double prev_b = colors[nl-1].b;

			double step_val = (val_min - prev_val_min)/(n_interp_col + 1);
			double step_r = ((double)r - prev_r)/(n_interp_col + 1);
			double step_g = ((double)g - prev_g)/(n_interp_col + 1);
			double step_b = ((double)b - prev_b)/(n_interp_col + 1);
			for( int i = 1; i <= n_interp_col; i++)
			{
				double val_min_interp = prev_val_min + i*step_val;
				int r_interp = (int)(prev_r + i*step_r);
				int g_interp = (int)(prev_g + i*step_g);
				int b_interp = (int)(prev_b + i*step_b);
				min_values.push_back(val_min_interp);
				HaColorMap::AddColor(r_interp,g_interp,b_interp);
			}
		}

		min_values.push_back(val_min);
		return HaColorMap::AddColor(r,g,b);
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in DValColorMap::AddColorAndMinVal() \n" );
		PrintLog("%s\n",ex.what());
		return -1;
	}
	return TRUE;
}

HaColor* DValColorMap::GetColorForVal(double val)
{
	if( !IsValid() )
	{
		PrintLog(" Error in DValColorMap::GetColorForVal() \n" );
		PrintLog( " Color map is not valid \n");
		return NULL;
	}

	int n = GetNColors();

	int i;
	for( i = n-1; i > 0; i-- )
	{
		if( val > min_values[i] ) return &colors[i];
	}
	return &colors[0];
}

int DValColorMap::SaveToTxtFile(const std::string& fname)
{
	try
	{
		char buf[256];
		if( !IsValid() ) throw std::runtime_error("Color Map is Not valid"); 
		std::ofstream ofs(fname.c_str());
		if(!ofs.is_open() ) throw std::runtime_error("Error to open text file to save Color Map"); 
		int n = GetNColors();
		int i;
		for( i = 0; i < n; i++ )
		{
			sprintf(buf, "%9.3e %3d %3d %3d", min_values[i], colors[i].r,colors[i].g,colors[i].b);
			ofs << buf << std::endl;
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in DValColorMap::SaveToTxtFile() \n");
		PrintLog("%s\n",ex.what());
	}
	return TRUE;
}

int DValColorMap::LoadFromTxtFile(const std::string& fname)
{
	try
	{
		std::ifstream ifs(fname.c_str());
		if(!ifs.is_open() ) throw std::runtime_error("Error to open text file to save Color Map"); 
		Clear();
		int r,g,b;
		double vmin;

		for(;;)
		{
			ifs >> vmin; if( ifs.fail() ) break;
			ifs >> r; if( ifs.fail() ) break;
			ifs >> g; if( ifs.fail() ) break;
			ifs >> b; if( ifs.fail() ) break;
			AddColorAndMinVal(r,g,b,vmin);
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in DValColorMap::LoadFromTxtFile() \n");
		PrintLog("%s\n",ex.what());
	}
	return TRUE;

}
