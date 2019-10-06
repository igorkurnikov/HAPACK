
include objlist

all:
	cd common && $(MAKE)
	cd object && $(MAKE)
	cd objfun && $(MAKE)
	cd ambfmt && $(MAKE)
	cd atmask && $(MAKE)
	cd capbox && $(MAKE)
	cd enefrc && $(MAKE)
	cd guilib && $(MAKE)
	cd format && $(MAKE)
	cd smarts && $(MAKE)
	cd tripos && $(MAKE)
	cd pdbent && $(MAKE)
	rm -f libmort.a
	ar cqs libmort.a ${OBJECTS}
clean:
	-(cd common && $(MAKE) clean)
	-(cd object && $(MAKE) clean)
	-(cd objfun && $(MAKE) clean)
	-(cd ambfmt && $(MAKE) clean)
	-(cd atmask && $(MAKE) clean)
	-(cd capbox && $(MAKE) clean)
	-(cd guilib && $(MAKE) clean)
	-(cd format && $(MAKE) clean)
	-(cd smarts && $(MAKE) clean)
	-(cd tripos && $(MAKE) clean)
	-(cd pdbent && $(MAKE) clean)
	rm -f libmort.a
	rm -f ${OBJECTS}
uninstall: clean
	-(cd common && $(MAKE) uninstall)
	-(cd object && $(MAKE) uninstall)
	-(cd objfun && $(MAKE) uninstall)
	-(cd ambfmt && $(MAKE) uninstall)
	-(cd atmask && $(MAKE) uninstall)
	-(cd capbox && $(MAKE) uninstall)
	-(cd guilib && $(MAKE) uninstall)
	-(cd format && $(MAKE) uninstall)
	-(cd smarts && $(MAKE) uninstall)
	-(cd tripos && $(MAKE) uninstall)
	-(cd pdbent && $(MAKE) uninstall)


