if INKSCAPE_EXT
INKSCAPE_DATADIR=$(datadir)@inkscape_ext_path@
EXTRA_DIST=extractpoly.inx extractpoly.py
install-data-local:
	$(MKDIR_P) $(INKSCAPE_DATADIR); \
	chmod u+w $(INKSCAPE_DATADIR); \
	$(INSTALL_DATA) $(top_srcdir)/extractpoly/extractpoly.inx $(INKSCAPE_DATADIR)/; \
	$(INSTALL_DATA) $(top_srcdir)/extractpoly/extractpoly.py $(INKSCAPE_DATADIR)/
uninstall-hook:
	-rm  $(INKSCAPE_DATADIR)/extractpoly.inx $(INKSCAPE_DATADIR)/extractpoly.py
endif
