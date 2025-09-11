########################
## Organizing man pages
########################

# nada

########################
## Class definitions
########################

#' Magma scenario class
#'
#' This entry documents the class definition used for [Magma]-type scenarios.
#' For details regarding the model, please refer to the [Magma model][Magma()] manual.
#'
#' @include class-Transferable.R
#' @slot growth_model *character*, selects the growth model, such as *exponential*
#'    or *logistic* growth.
#' @export
#' @seealso [Macrophyte-models]
setClass("Magma",
         contains=c("Transferable", "EffectScenario"),
         slots=c("growth_model"="character")
         )

# for backwards compatibility only
setClass("Myriophyllum", contains="Magma")
setClass("MyrioExp", contains="Myriophyllum")
setClass("MyrioExpScenario", contains="MyrioExp")
setClass("MyrioLog", contains="Myriophyllum")
setClass("MyrioLogScenario", contains="MyrioLog")


########################
## Constructor
########################

#' Magma model (Witt et al., submitted)
#'
#' The *Magma* model interprets the Tier 2C version of the *Lemna* model by
#' [Klein et al. (2021)][Lemna_SETAC()], as a generic macrophyte model.
#' It is mathematically equivalent to the Tier 2C version of the model
#' by Klein et al. (2021) with the recommended Tier 2C settings `k_photo_fixed=TRUE`
#' and `k_resp=0`.
#'
#' In particular, the growth model is a simple exponential growth model,
#' which is considered to be the typical situation for a laboratory macrophyte
#' study. Instead of frond numbers as for Lemna, the biomass is also returned as
#' total shoot length (*TSL*) in simulation results. Consequently, the model has
#' the additional parameter `r_DW_TSL` (dry weight per total shoot length ratio)
#' instead of `r_DW_FN` (dry weight per frond number ratio). A model variant
#' with an option for logistic growth is provided as well.
#'
#' @section State variables:
#' The model has two state variables:
#' - `BM`, Biomass (g dw)
#' - `M_int`, Mass of toxicant in plant population (ng)
#'
#' @section Model parameters:
#' The growth model can either simulate exponential growth (the default) or
#' logistic growth. For logistic growth, an additional parameter `D_L` describing
#' the limit density or carrying capacity needs to be provided.
#'
#' - Growth model
#'   - `mu_control`, Maximum photosynthesis rate (d-1), default: `0.47`
#'   - (optional) `D_L`, Limit density (g dw)
#'
#' - Concentration response (Toxicodynamics)
#'   - `EC50_int`, Internal concentration resulting in 50% effect (ug L-1)
#'   - `E_max`, Maximum inhibition (-), default: `1`
#'   - `b`, Slope parameter (-)
#'
#' - Internal concentration (Toxicokinetics)
#'   - `P`, Permeability (cm d-1)
#'   - `r_A_DW`, Area per dry-weight ratio (cm2 g-1), default: `1000`
#'   - `r_FW_DW`, Fresh weight per dry weight ratio (-), default: `16.7`
#'   - `r_FW_V`, Fresh weight density (g cm-3), default: `1`
#'   - `r_DW_TSL`, Dry weight per total shoot length ratio  (g dw cm-1)
#'   - `K_pw`, Partitioning coefficient plant:water (-), default: `1`
#'   - `k_met`, Metabolisation rate (d-1), default: `0`
#'
#' @section Environmental factors:
#'
#' None.
#'
#' @section Parameter boundaries:
#' Default values for parameter boundaries are set for all parameters by expert
#' judgement, for calibration purposes. Values can be modified using [set_bounds()].
#'
#' @section Simulation output:
#' Simulation results will contain the state variables biomass (`BM`) and
#' mass of internal toxicant (`M_int`).
#'
#' It is possible to amend the output of [simulate()] with additional model
#' quantities that are not state variables, for e.g. debugging purposes or to
#' analyze model behavior. To enable or disable additional outputs, use the
#' optional argument `nout` of [simulate()]. As an example, set `nout=2` to
#' enable reporting of total shoot length (`TSL`) and internal concentration
#' (`C_int`). Set `nout=0` to disable additional outputs. The default is `nout=1`.
#'
#' The available output levels are as follows:
#' - `nout` >= 1: `TSL`, total shoot length (cm)
#' - `nout` >= 2: `C_int`, internal concentration (ug L-1)
#' - `nout` >= 3: `f_photo`, photosynthesis dependency function (-)
#' - `nout` >= 4: `C_int_unb`, unbound internal concentration (ug L-1)
#' - `nout` >= 5: `C_ext`, external concentration (ug L-1)
#' - `nout` >= 6: `dBM`, biomass derivative (g dw d-1)
#' - `nout` >= 7: `dM_int`, mass of toxicant in plants derivative (ng d-1)
#'
#' @section Solver settings:
#' The arguments to ODE solver [deSolve::ode()] control how model equations
#' are numerically integrated. The settings influence stability of the numerical
#' integration scheme as well as numerical precision of model outputs. Generally, the
#' default settings as defined by *deSolve* are used, but all *deSolve* settings
#' can be modified in *cvasi* workflows by the user, if needed. Please refer
#' to e.g. [simulate()] on how to pass arguments to *deSolve* in *cvasi*
#' workflows.
#'
#' Some default settings of *deSolve* were adapted for this model by expert
#' judgement to enable precise, but also computationally efficient, simulations
#' for most model parameters. These settings can be modified by the user,
#' if needed:
#'
#' - `hmax = 0.1`<br>
#'    Maximum step length in time suitable for most simulations.
#'
#' @inheritSection Lemna_SETAC Effects
#' @inheritSection Transferable Biomass transfer
#' @param growth `character`, growth model to simulate: `"exp"` for exponential
#'    growth or `"log"` for logistic growth. Default is `"exp"`.
#' @references
#' Witt et al., submitted
#'
#' Klein J., Cedergreen N., Heine S., Reichenberger S., Rendal C.,
#' Schmitt W., Hommen U., 2021: *Refined description of the Lemna TKTD growth model
#' based on Schmitt et al. (2013) - equation system and default parameters*.
#' Report of the working group *Lemna* of the SETAC Europe Interest Group Effect
#' Modeling. Version 1.1, uploaded on 09 May 2022.
#' https://www.setac.org/group/effect-modeling.html
#'
#' @return an S4 object of type `Magma`
#' @seealso [Macrophyte-models], [Lemna-models], [Transferable], [Scenarios]
#' @family macrophyte models
#' @export
Magma <- function(growth=c("exp", "log")) {
  growth <- match.arg(growth)

  sc <- new("Magma",
      name=paste0("Magma-", growth),
      # `k_photo_max` only included for backwards compatibility of Myriophyllum-type scenarios
      param.req=c("mu_control", "E_max", "EC50_int", "b", "P", "r_A_DW",
                  "r_FW_DW", "r_FW_V", "r_DW_TSL", "K_pw", "k_met", "k_photo_max"),
      # default values as defined by Witt et al.
      param=list(mu_control=0.1, E_max=1, r_A_DW=1000, r_FW_DW=16.7, r_FW_V=1,
                 K_pw=1, k_met=0),
      # boundary presets defined by expert judgement
      param.bounds=list(mu_control=c(0, 1), E_max=c(0, 1), EC50_int=c(0, 1e6),
                        b=c(0.1, 20), P=c(0, 100)),
      endpoints=c("BM", "r"),
      control.req=TRUE,
      init=c(BM=0, M_int=0),
      transfer.comp.biomass="BM",
      transfer.comp.scaled="M_int",
      growth_model=growth
  ) %>% set_noexposure()

  # Adapt list of parameters for logistic growth
  if(growth == "log") {
    sc@param.req <- c(sc@param.req, "D_L")
    sc@param.bounds <- append(sc@param.bounds, list("D_L"=c(0, 1e6)))
  }

  sc
}

# constructors for backwards-compatibility
#' @export
#' @rdname Magma
Myrio <- function() {
  lifecycle::deprecate_soft(when="1.5.0", "Myrio()", "Magma()")
  Magma(growth="exp")
}

#' @export
#' @rdname Magma
Myrio_log <- function() {
  lifecycle::deprecate_soft(when="1.5.0", "Myrio_log()", "Magma(growth=\"log\")")
  Magma(growth="log")
}



########################
## Simulation
########################

# Numerically solve Magma scenarios
#
# @param scenario
# @param ... additional parameters passed on to [deSolve::ode()]
# @return data.frame
#' @include solver.R
#' @importFrom methods .hasSlot
solver_magma <- function(scenario, nout=1, method="lsoda", hmax=0.1, ...) {
  # for backwards compatibility with Myriophyllum-type scenarios
  if(!.hasSlot(scenario, "growth_model")) {
    # identify growth model to use
    gm <- "exp"
    if(is(scenario, "MyrioLog")) {
      gm <- "log"
    }
    # mock non-existing slot
    attr(scenario, "growth_model") <- gm

    # use `k_photo_max` as an alias of `mu_control`
    scenario@param["mu_control"] <- scenario@param["k_photo_max"]
  }
  # end of backwards compatibility

  params <- scenario@param

  # supplement parameter list with a placeholder value to ensure correct
  # parameter order
  if(scenario@growth_model == "exp") {
    params$growthno <- 1   # magic value: 1:=exp. growth, constant defined in model's C code
    params$D_L <- NA_real_ # magic value: dummy value w/o meaning, reqd to ensure parameter order
  } else {
    params$growthno <- 2 # magic value: 2:=log. growth, constant defined in model's C code
  }

  # make sure that parameters are present and in required order
  params_order <- c("mu_control", "growthno", "D_L", "E_max", "EC50_int", "b",
                    "P", "r_A_DW", "r_FW_DW", "r_FW_V", "r_DW_TSL", "K_pw",
                    "k_met")
  if(is.list(params)) {
    params <- unlist(params)
  }

  # check for missing parameters
  params_missing <- setdiff(params_order, names(params))
  if(length(params_missing)>0) {
    stop(paste("missing parameters:", paste(params_missing, collapse=",")))
  }

  # reorder parameters for deSolve
  params <- params[params_order]
  forcings <- list(scenario@exposure@series)

  # set names of additional output variables
  outnames <- c("TSL", "C_int", "f_photo", "C_int_unb", "C_ext", "dBM", "dM_int")

  ode(y=scenario@init, times=scenario@times, parms=params, dllname="cvasi",
        initfunc="magma_init", func="magma_func", initforc="magma_forc",
        forcings=forcings, nout=nout, method=method, hmax=hmax,
        outnames=outnames, ...)
}

#' @describeIn solver Numerically integrates `Magma` scenarios
setMethod("solver", "Magma", function(scenario, ...) solver_magma(scenario, ...))

########################
## Effects
########################

#' @include fx.R model-lemna_setac.R
#' @describeIn fx Effect at end of simulation of [Magma] scenarios
setMethod("fx", "Magma", function(scenario, ...) fx_lemna(scenario, ...))
