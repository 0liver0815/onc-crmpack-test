# This is to have along the study constant parameters settings of safety window length,
# whatever the cohort size is.
my_win_len <- SafetyWindowConst(
  patientGap = c(7, 5, 3),
  patientFollow = 7,
  patientFollowMin = 14
)
my_win_len
