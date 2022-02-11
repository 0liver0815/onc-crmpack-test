# Simulations-class.R line 36 add the following to prevent
# warning when package is loaded

# .GeneralSimulations <-
    # setClass(Class="GeneralSimulations",
             # representation(data="list",
                            # doses="numeric",
                            # seed="integer"),
             # prototype(data=
                           # list(Data(x=1:2,
                                     # y=0:1,
                                     # doseGrid=1:2,
                                     # #add to fix warning
                                     # ID=3:4, 
                                     # cohort=3:4),
                                # Data(x=3:4,
                                     # y=0:1,
                                     # doseGrid=3:4,
                                     # #add to fix warning
                                     # ID=3:4,
                                     # cohort=3:4)),
                       # doses=c(1, 2),
                       # seed=1L),
             # validity=
