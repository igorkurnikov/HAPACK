#ifndef PNPS_H
#define PNPS_H

#include <string>

/* Version number of package */
#define PNPSVERSION "8.3"

#define PNP_MAP_IO_STRING_LENGTH 512

class HaContext;

class PnpsObject
{
public:
	PnpsObject() { name = ""; }
	virtual ~PnpsObject() {}
	virtual int Clear() {return 0;}

	void SetName(const std::string& name_new) {
		name = name_new;
	}
	const char* GetName() const {
		return name.c_str();
	}
	const char* GetCStrName() const {
		return name.c_str();
	}
	std::string GetStdStrName() const {
		return name;
	}

protected:
	std::string name;
};
#endif
