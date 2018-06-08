/*!  \file ha_wx_aux_1.cpp

    Auxiliary Classes for HARLEM wxWidgets interface  

    \author Igor Kurnikov 
    \date 2010-

*/

#include <mpi.h>

#include "wx/wx.h"
#include "wx/valgen.h"

#include "haio.h"
#include "hatypes.h"
#include "ha_wx_aux_1.h"


#define HA_WX_AUX_1_CPP

wxDoubleValidator::wxDoubleValidator(double* fval_ptr_new, const wxString& format_new, double scale_new )
{
	fval_ptr = fval_ptr_new;
	format   = format_new;
	scale    = scale_new;
}

wxDoubleValidator::~wxDoubleValidator()
{

}
wxDoubleValidator::wxDoubleValidator(const wxDoubleValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool wxDoubleValidator::Copy(const wxDoubleValidator& val)
{
    wxValidator::Copy(val);

    fval_ptr = val.fval_ptr;
    format = val.format;
    scale = val.scale;

    return true;
}

wxObject* wxDoubleValidator::Clone() const
{
	wxDoubleValidator* ptr = new wxDoubleValidator(fval_ptr, format, scale);
	return ptr;
}

bool wxDoubleValidator::TransferToWindow()
{
	wxString str;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str.Printf(format.c_str(), (*fval_ptr) * scale );
	str.Trim(false);
	ctrl->SetValue(str);
	return true;
}

bool wxDoubleValidator::TransferFromWindow()
{
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	double tmp;
	bres = str.ToDouble(&tmp);
	if(!bres) return false;
	*fval_ptr = tmp/scale;
	return true;
}

bool wxDoubleValidator::Validate()
{
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	double tmp;
	bres = str.ToDouble(&tmp);
	if(!bres) return false;
	return true;
}


wxFloatValidator::wxFloatValidator(float* fval_ptr_new, const wxString& format_new)
{
	fval_ptr = fval_ptr_new;
	format   = format_new;
}

wxFloatValidator::~wxFloatValidator()
{

}
wxFloatValidator::wxFloatValidator(const wxFloatValidator& val)
    : wxValidator()
{
    Copy(val);
}
bool wxFloatValidator::Copy(const wxFloatValidator& val)
{
    wxValidator::Copy(val);

    fval_ptr = val.fval_ptr;
    format = val.format;

    return true;
}
wxObject* wxFloatValidator::Clone() const
{
	wxFloatValidator* ptr = new wxFloatValidator(fval_ptr, format);
	return ptr;
}

bool wxFloatValidator::TransferToWindow()
{
	wxString str;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str.Printf(format.c_str(),*fval_ptr);
	str.Trim(false);
	ctrl->SetValue(str);
	return true;
}

bool wxFloatValidator::TransferFromWindow()
{
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	double tmp;
	bres = str.ToDouble(&tmp);
	if(!bres) return false;
	*fval_ptr = tmp;
	return true;
}

bool wxFloatValidator::Validate()
{
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	double tmp;
	bres = str.ToDouble(&tmp);
	if(!bres) return false;
	return true;
}
//wxIntValidator
wxIntValidator::wxIntValidator(int* val_ptr_new, const wxString& format_new)
{
	//PrintLog("wxIntValidator::wxIntValidator\n");
	val_ptr = val_ptr_new;
	format   = format_new;
}

wxIntValidator::~wxIntValidator()
{

}
wxIntValidator::wxIntValidator(const wxIntValidator& val)
    : wxValidator()
{
    Copy(val);
}
bool wxIntValidator::Copy(const wxIntValidator& val)
{
    wxValidator::Copy(val);

    val_ptr = val.val_ptr;
    format = val.format;

    return true;
}
wxObject* wxIntValidator::Clone() const
{//PrintLog("wxIntValidator::Clone\n");
	wxIntValidator* ptr = new wxIntValidator(val_ptr, format);
	return ptr;
}

bool wxIntValidator::TransferToWindow()
{//PrintLog("wxIntValidator::TransferToWindow\n");
	wxString str;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str.Printf(format.c_str(),*val_ptr);
	str.Trim(false);
	ctrl->SetValue(str);
	return true;
}

bool wxIntValidator::TransferFromWindow()
{//PrintLog("wxIntValidator::TransferFromWindow\n");
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	long tmp;
	bres = str.ToLong(&tmp);
	if(!bres) return false;
	*val_ptr = (int)tmp;
	return true;
}

bool wxIntValidator::Validate()
{//PrintLog("wxIntValidator::Validate\n");
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	long tmp;
	bres = str.ToLong(&tmp);
	if(!bres) return false;
	return true;
}

//wxLongValidator
wxLongValidator::wxLongValidator(long* val_ptr_new, const wxString& format_new)
{
	//PrintLog("wxLongValidator::wxLongValidator\n");
	val_ptr = val_ptr_new;
	format   = format_new;
}

wxLongValidator::~wxLongValidator()
{

}
wxLongValidator::wxLongValidator(const wxLongValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool wxLongValidator::Copy(const wxLongValidator& val)
{
    wxValidator::Copy(val);

    val_ptr = val.val_ptr;
    format = val.format;

    return true;
}
wxObject* wxLongValidator::Clone() const
{//PrintLog("wxLongValidator::Clone\n");
	wxLongValidator* ptr = new wxLongValidator(val_ptr, format);
	return ptr;
}

bool wxLongValidator::TransferToWindow()
{//PrintLog("wxLongValidator::TransferToWindow\n");
	wxString str;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str.Printf(format.c_str(),*val_ptr);
	str.Trim(false);
	ctrl->SetValue(str);
	return true;
}

bool wxLongValidator::TransferFromWindow()
{//PrintLog("wxLongValidator::TransferFromWindow\n");
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	long tmp;
	bres = str.ToLong(&tmp);
	if(!bres) return false;
	*val_ptr = tmp;
	return true;
}

bool wxLongValidator::Validate()
{//PrintLog("wxLongValidator::Validate\n");
	wxString str;
	bool bres;
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	str = ctrl->GetValue();
	long tmp;
	bres = str.ToLong(&tmp);
	if(!bres) return false;
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// StdStringValidator
//
StdStringValidator::StdStringValidator(std::string* pstr_new)
{
	if(pstr_new == NULL)
	{
		PrintLog("Error in StdStringValidator::StdStringValidator() \n"); 
		PrintLog("Invalid String pointer \n");
	}
	pstr = pstr_new;
//	PrintLog(" pstr point to %s \n", (*pstr).c_str());
}

StdStringValidator::~StdStringValidator()
{

}
StdStringValidator::StdStringValidator(const StdStringValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool StdStringValidator::Copy(const StdStringValidator& val)
{
    wxValidator::Copy(val);

    pstr = val.pstr;

    return true;
}

wxObject* StdStringValidator::Clone() const
{
	StdStringValidator* ptr = new StdStringValidator(pstr);
	return ptr;
}

bool StdStringValidator::TransferToWindow()
{
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	ctrl->SetValue((*pstr).c_str());
	return true;
}

bool StdStringValidator::TransferFromWindow()
{
	wxTextCtrl* ctrl = (wxTextCtrl*) GetWindow();
	wxString str = ctrl->GetValue();
	*pstr = str.c_str();
	return true;
}

bool StdStringValidator::Validate()
{	
//	PrintLog("StdStringValidator::Validate() \n");
//	PrintLog(" pstr point to %s \n", (*pstr).c_str());
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// HaEnumValidator
//
HaEnumValidator::HaEnumValidator(HaEnum* p_enum_new)
{
	if(p_enum_new == NULL)
	{
		to_update_items = false;
		PrintLog("Error in HaEnumValidator::HaEnumValidator() \n"); 
		PrintLog("Invalid HaEnum pointer \n");
	}
	to_update_items = true;
	p_enum = p_enum_new;
}

HaEnumValidator::~HaEnumValidator()
{

}

HaEnumValidator::HaEnumValidator(const HaEnumValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool HaEnumValidator::Copy(const HaEnumValidator& val)
{
    wxValidator::Copy(val);
    p_enum = val.p_enum;
	to_update_items = val.to_update_items;
    return true;
}

wxObject* HaEnumValidator::Clone() const
{
	HaEnumValidator* ptr = new HaEnumValidator(p_enum);
	return ptr;
}

void HaEnumValidator::UpdateItems()
{
	wxControlWithItems* ctrl = (wxControlWithItems*) GetWindow();
	std::vector<std::string> str_arr = p_enum->GetActiveLabels();
	ctrl->Clear();
	int i;
	for( i=0; i < str_arr.size(); i++)
	{
		ctrl->Append(str_arr[i].c_str());
	}
	to_update_items = false;
}

bool HaEnumValidator::TransferToWindow()
{
	if( to_update_items ) UpdateItems();
	wxControlWithItems* ctrl = (wxControlWithItems*) GetWindow();
	ctrl->SetStringSelection(p_enum->label());
	return true;
}

bool HaEnumValidator::TransferFromWindow()
{
	wxControlWithItems* ctrl = (wxControlWithItems*) GetWindow();
	wxString str = ctrl->GetStringSelection();
	p_enum->SetWithLabel(str.c_str());
	return true;
}

bool HaEnumValidator::Validate()
{	
//	PrintLog("HaEnumValidator::Validate() \n");
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// HaEnumCheckBoxValidator
//
HaEnumCheckBoxValidator::HaEnumCheckBoxValidator(HaEnum* p_enum_new)
{
	if(p_enum_new == NULL)
	{
		PrintLog("Error in HaEnumCheckBoxValidator::HaEnumCheckBoxValidator() \n"); 
		PrintLog("Invalid HaEnum pointer \n");
	}
	p_enum = p_enum_new;
}

HaEnumCheckBoxValidator::~HaEnumCheckBoxValidator()
{

}

HaEnumCheckBoxValidator::HaEnumCheckBoxValidator(const HaEnumCheckBoxValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool HaEnumCheckBoxValidator::Copy(const HaEnumCheckBoxValidator& val)
{
    wxValidator::Copy(val);
    p_enum = val.p_enum;
    return true;
}

wxObject* HaEnumCheckBoxValidator::Clone() const
{
	HaEnumCheckBoxValidator* ptr = new HaEnumCheckBoxValidator(p_enum);
	return ptr;
}

bool HaEnumCheckBoxValidator::TransferToWindow()
{
	wxCheckBox* ctrl = (wxCheckBox*) GetWindow();
	ctrl->SetValue( p_enum->value() != 0 );
	return true;
}

bool HaEnumCheckBoxValidator::TransferFromWindow()
{
	wxCheckBox* ctrl = (wxCheckBox*) GetWindow();
	int val = ctrl->GetValue();
	p_enum->SetWithValue(val);
	return true;
}

bool HaEnumCheckBoxValidator::Validate()
{	
//	PrintLog("HaEnumValidator::Validate() \n");
	return true;
}

//////////////////////////////////////////////////////////////////////////////
// IntCheckBoxValidator
//
IntCheckBoxValidator::IntCheckBoxValidator(int* p_int_new)
{
	if(p_int_new == NULL)
	{
		PrintLog("Error in IntCheckBoxValidator::IntCheckBoxValidator() \n"); 
		PrintLog("Invalid int pointer \n");
	}
	p_int = p_int_new;
}

IntCheckBoxValidator::~IntCheckBoxValidator()
{

}

IntCheckBoxValidator::IntCheckBoxValidator(const IntCheckBoxValidator& val)
    : wxValidator()
{
    Copy(val);
}

bool IntCheckBoxValidator::Copy(const IntCheckBoxValidator& val)
{
    wxValidator::Copy(val);
    p_int = val.p_int;
    return true;
}

wxObject* IntCheckBoxValidator::Clone() const
{
	IntCheckBoxValidator* ptr = new IntCheckBoxValidator(p_int);
	return ptr;
}

bool IntCheckBoxValidator::TransferToWindow()
{
	wxCheckBox* ctrl = (wxCheckBox*) GetWindow();
	ctrl->SetValue( *p_int != 0 );
	return true;
}

bool IntCheckBoxValidator::TransferFromWindow()
{
	wxCheckBox* ctrl = (wxCheckBox*) GetWindow();
	(*p_int) = ctrl->GetValue() ? 1 : 0;
	return true;
}

bool IntCheckBoxValidator::Validate()
{	
//	PrintLog("IntCheckBoxValidator::Validate() \n");
	return true;
}

void DDX_Choice_HaEnum(wxWindow* win, int id, HaEnum* p_enum)
{
	wxChoice* choice_ctrl = (wxChoice*) win->FindWindow(id);
	choice_ctrl->SetValidator( HaEnumValidator(p_enum) );
}

void DDX_Text_double(wxWindow* win, int id, double& dbl, const wxString& format)
{
	double* ptr_dbl = &dbl;
	wxTextCtrl* edit_ctrl;
	edit_ctrl = (wxTextCtrl*) win->FindWindow(id);
	edit_ctrl->SetValidator( wxDoubleValidator(ptr_dbl, format) );
}

void DDX_Text_int(wxWindow* win, int id, int& int_num)
{
	int* ptr_int = &int_num;
	wxTextCtrl* edit_ctrl;
	edit_ctrl = (wxTextCtrl*) win->FindWindow(id);
	edit_ctrl->SetValidator( wxGenericValidator(ptr_int) );
}

void DDX_Text_long(wxWindow* win, int id, long& long_num)
{
	long* ptr_long = &long_num;
	wxTextCtrl* edit_ctrl;
	edit_ctrl = (wxTextCtrl*) win->FindWindow(id);
	edit_ctrl->SetValidator( wxLongValidator(ptr_long,"%d") );
}

void AddWindowToArray(wxWindow* parent, int id, std::vector<wxWindow*>& win_arr)
{
	wxWindow* ctrl_win = parent->FindWindow(id);
	if( ctrl_win != NULL ) 
	{
		win_arr.push_back(ctrl_win);
	}
	else
	{
		PrintLog("Error in AddCtrlToArray() \n");
		PrintLog("Window with ID = %d  not found \n",id);
	}
}