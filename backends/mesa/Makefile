include $(MESA_DIR)/utils/makefile_header
.PHONY: all
all:
	$(FC) $(FCopenmp) -O2 -fPIC -ffree-line-length-none -I$(MESA_DIR)/include -J$(DEST) -dynamiclib $(SOURCE) -o $(DEST)/libqstars_mesa.dylib -L$(MESA_DIR)/lib -latm -lturb -lkap -leos -lchem $(LOAD_MESA_NUMERICS)
	$(FC) $(FCopenmp) -O2 -ffree-line-length-none -I$(DEST) -I$(MESA_DIR)/include $(dir $(SOURCE))check_microphysics.f90 $(DEST)/libqstars_mesa.dylib -o $(DEST)/check_microphysics
