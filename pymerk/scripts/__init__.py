from pydantic import BaseModel, Field


# from https://github.com/grimme-lab/CENSO/blob/main/example.censo2rc
class General(BaseModel):
    temperature: float = Field(298.15, description='Default temperature in Kelvin used for xtb and COSMOtherm')
    evaluate_rrho: bool = Field(True, description='Calculate thermal energy contributions using xtb')
    sm_rrho: str = Field('gbsa', description='Solvation model used for rrho corrections')
    imagthr: float = Field(-100.0, description='Threshold for accepting imaginary frequencies (cm^-1)')
    sthr: float = Field(50.0, description='Wave number threshold for switching in rrho approximation (cm^-1)')
    solvent: str = Field('h2o', description='Default solvent identifier (used by solvation models)')
    gas_phase: bool = Field(False, description='If True, calculations treat system as gas-phase')
    # copy_mo: bool = Field(True, description='Copy molecular orbitals between steps when supported')
    # balance: bool = Field(True, description='Attempt to maximize core utilization')
    # ignore_failed: bool = Field(True, description='Continue workflow even if intermediate jobs fail')


class Prescreening(BaseModel):
    prog: str = Field('vlx', description='Program used for prescreening jobs')
    func: str = Field('pbe0', description='Functional used for prescreening')
    basis: str = Field('def2-svp', description='Basis set for prescreening')
    gfnv: str = Field('gfn2', description='xTB variant for gsolv contributions')
    threshold: float = Field(4.0, description='Energy threshold (kcal/mol) to keep candidates')
    # template: bool = Field(False, description='Tries to insert template file')


class Screening(BaseModel):
    prog: str = Field('vlx', description='Program used for screening')
    func: str = Field('rcam-b3lyp', description='Functional used for screening')
    basis: str = Field('def2-tzvpd', description='Basis set used during screening')
    sm: str = Field('smd', description='Solvation model for screening')
    gfnv: str = Field('gfn2', description='xTB variant')
    threshold: float = Field(3.5, description='Energy cutoff for keeping candidates (kcal/mol)')
    gsolv_included: bool = Field(
        False, description='Whether solvation free energy should be included in energies or calculated separately')
    # template: bool = Field(False, description='Tries to insert template file')


class Optimization(BaseModel):
    prog: str = Field('vlx', description='Program used for optimizations')
    func: str = Field('rcam-b3lyp', description='Functional used for optimizations')
    basis: str = Field('def2-dzvp', description='Basis set for optimizations')
    sm: str = Field('cpcm', description='Solvation model')
    gfnv: str = Field('gfn2', description='xTB variant for any semiempirical steps')
    optcycles: int = Field(8, description='Number of optimization macrocycles to attempt')
    maxcyc: int = Field(200, description='Maximum optimization microcycles')
    optlevel: str = Field('normal', description='Optimization thoroughness (e.g., loose, normal, tight)')
    threshold: float = Field(3.0, description='Energy threshold (kcal/mol)')
    gradthr: float = Field(
        0.01, description='Gradient threshold (a.u.) below which the energy threshold will be applied')
    # hlow: float = Field(0.01, description='Low-level threshold parameter')
    macrocycles: bool = Field(True, description='Allow macrocycle protocol')
    # constrain: bool = Field(False, description='Apply geometry constraints if True')
    # xtb_opt: bool = Field(True, description='Whether to use ANCOPT as driver')
    # template: bool = Field(False, description='Tries to insert template file')


class Refinement(BaseModel):
    prog: str = Field('vlx', description='Program used for refinement')
    func: str = Field('wb97m-d4', description='Functional used for refinement')
    basis: str = Field('def2-tzvp', description='Basis set for refinement')
    sm: str = Field('cpcm', description='Solvation model for refinement')
    gfnv: str = Field('gfn2', description='xTB variant')
    threshold: float = Field(0.95, description='Boltzmann population threshold')
    # gsolv_included: bool = Field(
    #     False, description='Whether solvation free energy should be included in energies or calculated separately')
    # template: bool = Field(False, description='Tries to insert template file')


class Paths(BaseModel):
    xtb: str = Field('', description='xtb binary path')
    vlx: str = Field('', description='VeloxChem binary path')


class Config(BaseModel):
    general: General = Field(default_factory=General)
    prescreening: Prescreening = Field(default_factory=Prescreening)
    screening: Screening = Field(default_factory=Screening)
    optimization: Optimization = Field(default_factory=Optimization)
    refinement: Refinement = Field(default_factory=Refinement)
    paths: Paths = Field(default_factory=Paths)
