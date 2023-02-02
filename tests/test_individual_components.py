# DefaultPolicyMakerModel
# Invalid filepath
# Malformated file


# DefaultSeasonalEffectsModel
# seasonal_multiplier_by_region_by_month is null (this is valid)
# seasonal_multiplier_by_region_by_month points to an non-existant file (this is invalid)
# Mismatch between regions in file and regions defined elsewhere


# DefaultVaccinationModel
# Why are `number_of_strains` and `number_of_rho_immunity_outcomes` required, given the values can be calculated from values in the health model?
# Probably better to make this a ValueError rather than an AssertionError.
