import argparse
import pathlib
import shutil
import tomllib

import pymerk
from pymerk.driver import XtbDriver, VlxDriver
from pymerk.ensemble import Ensemble
from pymerk.scripts import Config
from pymerk.scripts.filter import EnergyFilter


class BaseWorkflow:
    def __init__(self, workdir: str | pathlib.Path):
        self.workdir = pathlib.Path(workdir)

    def filter(self, ensemble: Ensemble) -> Ensemble:
        return ensemble


def get_xtb_driver(config: Config, workdir: pathlib.Path, **kwargs) -> XtbDriver:
    if config.paths.xtb == '':
        raise RuntimeError('No `xtb` path specified')

    opt_dict = dict(
        imagthr=config.general.imagthr,
        sthr=config.general.sthr,
        scale=config.general.scale
    )

    if not config.general.gas_phase:
        opt_dict.update(solvatation_model=config.general.sm_rrho, solvent=config.general.solvent)

    return XtbDriver(
        workdir, config.paths.xtb,
        **kwargs
    )


def get_vlx_driver(config: Config, workdir: pathlib.Path, **kwargs) -> VlxDriver:
    if config.paths.vlx == '':
        raise RuntimeError('No `vlx` path specified')

    return VlxDriver(workdir, config.paths.vlx, **kwargs)


GET_DRIVER = {
    'xtb': get_xtb_driver,
    'vlx': get_vlx_driver
}


def hp(title: str):
    print()
    print('*' * (len(title) + 4))
    print('* {} *'.format(title))
    print('*' * (len(title) + 4))
    print()


class DefaultWorkflow(BaseWorkflow):
    """A censo-like workflow, with prescreening → screening → optimisation → refinement
    """

    AU_TO_KCAL = 6.275030e2

    def __init__(self, workdir: str | pathlib.Path, config: Config):
        super().__init__(workdir)
        self.config = config

    def filter(self, ensemble: Ensemble) -> Ensemble:
        print('* Input: {} conformer(s)'.format(len(ensemble)))

        hp('Prescreening')  # TODO: missing gsolv, technically at the gbsa level

        workdir = self.workdir / '1_prescreening'
        workdir.mkdir(exist_ok=True)
        output_file = self.workdir / '1_preescreening.log'

        with output_file.open('w') as f:
            new_ensemble = EnergyFilter(
                GET_DRIVER[self.config.prescreening.prog](
                    self.config, workdir,
                    method=self.config.prescreening.func,
                    basis=self.config.prescreening.basis),
                self.config.prescreening.threshold / self.AU_TO_KCAL
            ).filter(ensemble, f)

        shutil.rmtree(workdir)

        return new_ensemble


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=pathlib.Path, help='Input file')
    parser.add_argument('conformers', type=pathlib.Path, help='XYZ file with the conformers')
    parser.add_argument(
        '-o', '--output', type=pathlib.Path, default='final_ensemble.xyz', help='final set of conformers')

    parser.add_argument('-c', '--charge', type=int, default=0, help='Charge of the conformers')
    parser.add_argument('-m', '--multiplicity', type=int, default=1, help='Multiplicity of the conformers')

    args = parser.parse_args()

    print('* This is pymerk v{}'.format(pymerk.__version__))

    # load config
    config = Config()

    if args.input:
        with args.input.open('rb') as f:
            config = Config(**tomllib.load(f))

    # load ensemble
    with args.conformers.open('r') as f:
        conformers = Ensemble.from_multi_xyz(f, args.charge, args.multiplicity)

    # sort
    new_conformers = DefaultWorkflow(pathlib.Path.cwd(), config).filter(conformers)

    # save
    with args.output.open('w') as f:
        new_conformers.as_multi_xyz(f)


if __name__ == '__main__':
    main()
