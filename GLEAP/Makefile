include ../config.h

install:  sleap

sleap:
	@echo You must compile sleap/gleap with GNU make

clean:
	cd freelib && $(MAKE) clean
	cd mortsrc && $(MAKE) clean
	cd plugins && $(MAKE) clean
	cd leapsrc && $(MAKE) clean
	cd mortsrc/mortgl && $(MAKE) clean
	cd leapsrc/widget && $(MAKE) clean

gleap: install
	cd mortsrc/mortgl && $(MAKE)
	cd leapsrc/widget && $(MAKE) ; cp gleap $(BINDIR)

uninstall:
	-(cd freelib && $(MAKE) uninstall)
	-(cd mortsrc && $(MAKE) uninstall)
	-(cd plugins && $(MAKE) uninstall)
	-(cd leapsrc && $(MAKE) uninstall)
	-(cd mortsrc/mortgl && $(MAKE) uninstall)
	-(cd leapsrc/widget && $(MAKE) uninstall)
	-rm -f $(BINDIR)/gleap
