mcMd_modules_micelle_analyzers_= \
    mcMd/modules/micelle/analyzers/ClusterHistogramInd.cpp \
    mcMd/modules/micelle/analyzers/ClusterHistogramSimple.cpp \
    mcMd/modules/micelle/analyzers/ClusterIdentifierSG.cpp \
    mcMd/modules/micelle/analyzers/InterfacialLoading.cpp \
    mcMd/modules/micelle/analyzers/MicelleFlux.cpp \
    mcMd/modules/micelle/analyzers/MicelleFluxDroplet.cpp \
    mcMd/modules/micelle/analyzers/RadialComposition.cpp

mcMd_modules_micelle_analyzers_SRCS=\
  $(addprefix $(SRC_DIR)/, $(mcMd_modules_micelle_analyzers_))

mcMd_modules_micelle_analyzers_OBJS=\
  $(addprefix $(BLD_DIR)/, $(mcMd_modules_micelle_analyzers_:.cpp=.o))

