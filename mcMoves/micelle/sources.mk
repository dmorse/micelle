mcMd_modules_micelle_mcMoves_micelle_=\
    mcMd/mcMoves/micelle/AggregatorMove.cpp \
    mcMd/mcMoves/micelle/SingleMicelleHybridMove.cpp \
    mcMd/mcMoves/micelle/SingleMicelleUmbrellaSamplingMove.cpp \
    mcMd/mcMoves/micelle/ClusterIdentifierMC.cpp \
    mcMd/mcMoves/micelle/RgUmbrellaMove.cpp

mcMd_modules_micelle_mcMoves_micelle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_modules_micelle_mcMoves_micelle_))
mcMd_modules_micelle_mcMoves_micelle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_modules_micelle_mcMoves_micelle_:.cpp=.o))

