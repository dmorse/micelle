mcMd_modules_micelle_mcMoves_micelle_=\
    mcMd/modules/micelle/mcMoves/micelle/SingleMicelleHybridMove.cpp \
    mcMd/modules/micelle/mcMoves/micelle/SingleMicelleUmbrellaSamplingMove.cpp \
    mcMd/modules/micelle/mcMoves/micelle/ClusterIdentifierMC.cpp \
    mcMd/modules/micelle/mcMoves/micelle/RgUmbrellaMove.cpp

mcMd_modules_micelle_mcMoves_micelle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_modules_micelle_mcMoves_micelle_))
mcMd_modules_micelle_mcMoves_micelle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_modules_micelle_mcMoves_micelle_:.cpp=.o))

