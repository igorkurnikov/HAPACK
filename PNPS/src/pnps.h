#ifndef PNPSOBJECT_H
#define PNPSOBJECT_H

#include <string>

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
