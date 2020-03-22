noinst_LIBRARIES = libipack.a

AM_CXXFLAGS = -O3 -Dlinux -DGNU -DGCC_NEW -DTIMER -fPIC 

AM_CPPFLAGS = -I$(top_srcdir)/WOLFGANG/generic \
 -I$(top_srcdir)/WOLFGANG/IPACK \
 -I$(top_srcdir)/WOLFGANG/DO_LIB \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/include \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/template \
  $(MPI_INCLUDE)
libipack_a_SOURCES= ipack_main.cpp \
      function.cpp  integ_array.cpp coef_array.cpp rec.cpp basis.cpp f.cpp \
      quad.cpp contract.cpp vrr.cpp twoel0.cpp hrr1.cpp kinetic.cpp overlap.cpp nuclear.cpp \
      coef_set.cpp hrr2.cpp  input.cpp operators.cpp  symdesig.cpp\
      hrr_set.cpp storage.cpp output1.cpp output2.cpp ../generic/parallel.cpp print.cpp output3.cpp \
      sphere.cpp symop.cpp symlogic.cpp twoel.cpp rec_stor.cpp one_el.cpp \
      four_stor.cpp const.cpp orbinfo.cpp types.cpp nuc_rep.cpp ../generic/vault.cpp \
      element.cpp e_field.cpp moment.cpp angular.cpp spin_orb.cpp alloc.cpp fix.cpp \
      integ_file.cpp nuc_mom.cpp nuc_e_field.cpp qc_utilities.cpp
