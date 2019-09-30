include ../config.h

# BUILD_SLEAP is defined to either install_sleap or skip_sleap by configure_at
install:  $(BUILD_SLEAP)

install_sleap:  sleap

sleap:
	cd freelib && $(MAKE)
	cd mortsrc && $(MAKE)
	cd plugins && $(MAKE)
	cd leapsrc && $(MAKE) install

skip_sleap:
	@echo "NOTE: skipping the building of sleap because the compiler version is"
	@echo "      too old or -nosleap was specified as an option to configure."

clean:
	-(cd freelib && $(MAKE) clean)
	-(cd mortsrc && $(MAKE) clean)
	-(cd plugins && $(MAKE) clean)
	-(cd leapsrc && $(MAKE) clean)
	-(cd mortsrc/mortgl && $(MAKE) clean)
	-(cd leapsrc/widget && $(MAKE) clean)

gleap: install
	cd mortsrc/mortgl && $(MAKE)
	cd leapsrc/widget && $(MAKE) ; cp gleap $(BINDIR)

uninstall: clean
	-(cd freelib && $(MAKE) uninstall)
	-(cd mortsrc && $(MAKE) uninstall)
	-(cd plugins && $(MAKE) uninstall)
	-(cd leapsrc && $(MAKE) uninstall)
	-(cd mortsrc/mortgl && $(MAKE) uninstall)
	-(cd leapsrc/widget && $(MAKE) uninstall)
	-rm -f $(BINDIR)/gleap
