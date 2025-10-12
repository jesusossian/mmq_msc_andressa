NAME          SPMETALS
ROWS
 N  VALUE
 N  CHANGOBJ
 E  YIELD
 L  FE
 L  MN
 L  CU
 L  MG
 G  AL
 G  SI
COLUMNS
    BIN1      VALUE             .03    YIELD         1.00
    BIN1      FE                .15    MN             .02
    BIN1      CU                .03    MG             .02
    BIN1      AL                .30    SI              .5
    BIN1      CHANGOBJ          .09
    BIN2      VALUE             .08    YIELD         1.00
    BIN2      FE                .04    MN             .04
    BIN2      CU                .05    MG             .03
    BIN2      AL                .75    SI             .06
    BIN3      VALUE             .17    YIELD         1.00
    BIN3      FE                .02    MN             .01
    BIN3      CU                .08    AL             .80
    BIN3      SI                .08
    BIN3      CHANGOBJ          .25
    BIN4      VALUE             .12    YIELD         1.00
    BIN4      FE                .04    MN             .02
    BIN4      CU                .02    AL             .75
    BIN4      SI                .12
    BIN5      VALUE             .15    YIELD         1.00
    BIN5      FE                .02    MN             .02
    BIN5      CU                .06    MG             .01
    BIN5      AL                .80    SI             .02
    ALUM      VALUE             .21    YIELD         1.00
    ALUM      FE                .01    CU             .01
    ALUM      AL                .97    SI             .01
    SILICON   VALUE             .38    YIELD         1.00
    SILICON   FE                .03    SI             .97
RHS
    RHS       YIELD           2000.    FE             60.
    RHS       CU               100.    MN             40.
    RHS       MG                30.    AL            1500.
    RHS       SI               300.
    CHANGRHS  SI               100.
    CHANGRHS  CU                60.
RANGES
    RNG       SI                50.
    CHANGRNG  SI                30.
    CHANGRNG  CU                40.
BOUNDS
 UP BNN       BIN1            200.
 UP BNN       BIN2            750.
 LO BNN       BIN3            400.
 UP BNN       BIN3            800.
 LO BNN       BIN4            100.
 UP BNN       BIN4            700.
 UP BNN       BIN5            1500.
 UP CHANGBND  BIN4            200.
 UP CHANGBND  BIN5            1800.
ENDATA

