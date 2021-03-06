
# In this example we define a rule for dose increments which would allow:
# Maximum skip one dose level, that is 2 dose levels higher than the last dose
# given.
# Maximum increment is explicitly defined as:
increments <- IncrementsNumDoseLevels(maxLevels = 2, basisLevel = "last")
# Since the default method is based on the last dose given,
# maximum increment can also be defined as:
increments <- IncrementsNumDoseLevels(maxLevels = 2)
