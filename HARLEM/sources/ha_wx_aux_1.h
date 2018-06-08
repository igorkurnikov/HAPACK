/*!  \file ha_wx_aux_1.h

    Auxiliary Classes for HARLEM wxWidgets interface  

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef HA_WX_AUX_1_H
#define HA_WX_AUX_1_H

#include <vector>

#include "wx/validate.h"


class wxDoubleValidator: public wxValidator
//! Validator for double numbers
{
public: 
	wxDoubleValidator(double* fval_ptr_new, const wxString& format_new, double scale_new = 1.0);
	virtual ~wxDoubleValidator();
    wxDoubleValidator(const wxDoubleValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
    bool Copy(const wxDoubleValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4

	double*  fval_ptr;
	wxString format;
	double   scale; // factor to multiply the value when displayed in the window( from Bohr to Ang for example)

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
private:
// Cannot use
//  DECLARE_NO_COPY_CLASS(wxTextValidator)
// because copy constructor is explicitly declared above;
// but no copy assignment operator is defined, so declare
// it private to prevent the compiler from defining it:
//    wxDoubleValidator& operator=(const wxDoubleValidator&);
};

class wxFloatValidator: public wxValidator
//! Validator for float numbers
{
public: 
	wxFloatValidator(float* fval_ptr_new, const wxString& format_new);
	virtual ~wxFloatValidator();
  wxFloatValidator(const wxFloatValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
  bool Copy(const wxFloatValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
	float*  fval_ptr;
	wxString format;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

class wxIntValidator: public wxValidator
//! Validator for int numbers
{
public: 
	wxIntValidator(int* val_ptr_new, const wxString& format_new);
	virtual ~wxIntValidator();
  wxIntValidator(const wxIntValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
  bool Copy(const wxIntValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
	int*  val_ptr;
	wxString format;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

class wxLongValidator: public wxValidator
//! Validator for int numbers
{
public: 
	wxLongValidator(long* val_ptr_new, const wxString& format_new);
	wxLongValidator(const wxLongValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
	virtual ~wxLongValidator();
 
    bool Copy(const wxLongValidator& val);//mikola 09-13-2007 compatibility with gcc >=3.4
	long*  val_ptr;
	wxString format;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

class StdStringValidator: public wxValidator
//! Validator for std::string
{
public: 
	StdStringValidator(std::string* pstr_new);
	virtual ~StdStringValidator();
    StdStringValidator(const StdStringValidator& val);  //mikola 09-13-2007 compatibility with gcc >=3.4
    bool Copy(const StdStringValidator& val);           //mikola 09-13-2007 compatibility with gcc >=3.4
	std::string*  pstr;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

class HaEnum;

class HaEnumValidator: public wxValidator
//! Validator for HaEnum object and wxChoice widget 
{
public: 
	HaEnumValidator(HaEnum* p_enum);
	virtual ~HaEnumValidator();
    HaEnumValidator(const HaEnumValidator& val);  
    bool Copy(const HaEnumValidator& val);          
	HaEnum*  p_enum;

	void UpdateItems();      //!< Update list of items in control associated with Validator
	void SetToUpdateItems(bool set_par = true ) { to_update_items = set_par; }
	bool ToUpdateItems() const { return to_update_items; }
	bool to_update_items; //!< flag to update list of items in control associated with Validator

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

class HaEnumCheckBoxValidator: public wxValidator
//! Validator for HaEnum object (should have only 0 and 1 enum values ) and wxCheckBox widget 
{
public: 
	HaEnumCheckBoxValidator(HaEnum* p_enum);
	virtual ~HaEnumCheckBoxValidator();
    HaEnumCheckBoxValidator(const HaEnumCheckBoxValidator& val);  
    bool Copy(const HaEnumCheckBoxValidator& val);          
	HaEnum*  p_enum;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};


class IntCheckBoxValidator: public wxValidator
//! Validator for int(TRUE or FALSE) and wxCheckBox widget 
{
public: 
	IntCheckBoxValidator(int* p_int);
	virtual ~IntCheckBoxValidator();
    IntCheckBoxValidator(const IntCheckBoxValidator& val);  
    bool Copy(const IntCheckBoxValidator& val);          
	int*  p_int;

	virtual wxObject* Clone() const;
	virtual bool TransferToWindow();
	virtual bool TransferFromWindow();
	virtual bool Validate();
};

#if !defined(HA_WX_AUX_1_CPP)
  extern void DDX_Text_double(wxWindow* win, int id, double& dbl, const wxString& format);
  extern void DDX_Text_int(wxWindow* win,  int id, int& int_num);
  extern void DDX_Text_long(wxWindow* win, int id, long& long_num);
  extern void DDX_Choice_HaEnum(wxWindow* win, int id, HaEnum* p_enum);
  extern void AddWindowToArray(wxWindow* parent, int id, std::vector<wxWindow*>& win_arr);
#else
  // void InitWxChoiceHaEnum(wxWindow* win, HaEnum* p_enum);
  void DDX_Text_double(wxWindow* win, int id, double& dbl, const wxString& format);
  void DDX_Text_int(wxWindow* win, int id, int& int_num);
  void DDX_Text_long(wxWindow* win, int id, long& long_num);
#endif

#endif // End of #define HA_WX_AUX_1_H