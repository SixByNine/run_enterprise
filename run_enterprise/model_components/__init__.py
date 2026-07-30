from . import red_model, white_model, dm_model, timing_model, glitch_model, planet_model, solarwind_model, \
    chromatic_model, nudot_switch_model
from . import chol_red_model, chol_white_model

all = [red_model, white_model, dm_model, timing_model, glitch_model, planet_model, solarwind_model, chromatic_model,nudot_switch_model]

try:
    from . import swgp
except ModuleNotFoundError as exc:
    # Keep the package importable when optional enterprise_extensions is absent.
    if not exc.name or not exc.name.startswith("enterprise_extensions"):
        raise
else:
    all.append(swgp)

cov_models = [chol_red_model, chol_white_model]
res_models = [timing_model, glitch_model]
