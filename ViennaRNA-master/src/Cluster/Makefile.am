AUTOMAKE_OPTIONS = subdir-objects

include $(top_srcdir)/silent_rules.mk

bin_PROGRAMS = AnalyseDists AnalyseSeqs
man_MANS = AnalyseDists.1 AnalyseSeqs.1

AM_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_builddir)
AM_LDFLAGS = $(RNA_LDFLAGS) $(PTHREAD_LIBS)

if VRNA_AM_SWITCH_STATIC_EXECUTABLES
AM_LDFLAGS += -all-static
endif

LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la

# Link against stdc++ if we use SVM
if VRNA_AM_SWITCH_SVM
LDADD += $(SVM_LIBS)
endif

AnalyseDists_SOURCES = split.c \
          distance_matrix.c \
          cluster.c \
          treeplot.c \
          AD_main.c

AnalyseSeqs_SOURCES = statgeom.c \
          distance_matrix.c \
          cluster.c \
          AS_main.c \
          split.c \
          treeplot.c \
          PS3D.c

noinst_HEADERS = cluster.h split.h treeplot.h distance_matrix.h statgeom.h \
                PS3D.h StrEdit_CostMatrix.h  

EXTRA_DIST = $(man_MANS)
