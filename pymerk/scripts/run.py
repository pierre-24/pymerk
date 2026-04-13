import argparse
import pathlib
import shutil
import tomllib
from typing import Any

import pymerk
from pymerk.driver import XtbDriver, VlxDriver, BaseDriver
from pymerk.ensemble import Ensemble
from pymerk.scripts import Config
from pymerk.scripts.filter import (
    EnergyFilter, EnergyWithXtbGsolvFilter, GibbsFreeEnergyWithXtbFilter, GibbsFreeEnergyFilter, BaseFilter)

# Constants
AU_TO_KCAL = 6.275030e2


class BaseWorkflow:
    def __init__(self, workdir: str | pathlib.Path, config: Any):
        self.workdir = pathlib.Path(workdir)
        self.config = config

    @staticmethod
    def _print_header(title: str):
        border = '*' * (len(title) + 4)
        print(f'\n{border}\n* {title} *\n{border}\n')

    def _get_driver(self, prog: str, stage_dir: pathlib.Path, **overrides) -> BaseDriver:
        """Centralized driver factory."""
        paths = self.config.paths

        if prog == 'xtb':
            if not paths.xtb:
                raise RuntimeError('xtb path not configured.')

            # Default XTB options from general config
            opts = {
                'imagthr': self.config.general.imagthr,
                'sthr': self.config.general.sthr,
            }
            if not self.config.general.gas_phase:
                opts.update({
                    'solvatation_model': self.config.general.sm_rrho,
                    'solvent': self.config.general.solvent
                })
            # Apply method/version overrides from stage or aux calls
            opts.update(overrides)
            return XtbDriver(stage_dir, paths.xtb, **opts)

        if prog == 'vlx':
            if not paths.vlx:
                raise RuntimeError('VeloxChem path not configured.')

            return VlxDriver(stage_dir, paths.vlx, **overrides)

        raise ValueError(f'Unsupported driver: {prog}')

    def _get_filter(self, stage_cfg: Any, stage_dir: pathlib.Path, main_driver: BaseDriver, f_type: str) -> BaseFilter:
        """Filter factory that decides which logic to apply based on config."""

        threshold = stage_cfg.threshold / AU_TO_KCAL

        if f_type == 'energy':
            return EnergyFilter(main_driver, threshold)

        if f_type == 'gsolv':
            # Use _get_driver to create the auxiliary XTB driver
            xtb_aux = self._get_driver('xtb', stage_dir, version=stage_cfg.gfnv)
            return EnergyWithXtbGsolvFilter(main_driver, xtb_aux, threshold)

        if f_type == 'gibbs':
            temp = self.config.general.temperature
            return GibbsFreeEnergyFilter(main_driver, threshold, T=temp)

        if f_type == 'gibbs_xtb':
            xtb_aux = self._get_driver('xtb', stage_dir, version=stage_cfg.gfnv)
            temp = self.config.general.temperature
            return GibbsFreeEnergyWithXtbFilter(main_driver, xtb_aux, threshold, T=temp)

        raise ValueError(f'Unknown filter type: {f_type}')

    def _run_stage(self, name: str, stage_cfg: Any, ensemble: Ensemble, f_type: str = float) -> Ensemble:
        """Generic stage runner: handles IO, driver init, and filtering."""
        BaseWorkflow._print_header(f'Stage: {name}')

        stage_dir = self.workdir / name
        log_file = self.workdir / f'{name}.log'
        stage_dir.mkdir(exist_ok=True)

        try:
            with log_file.open('w') as f:
                # 1. Create Main Driver
                opts = {}

                if hasattr(stage_cfg, 'sm') and stage_cfg.sm != '':  # add solvent if any
                    opts.update(solvatation_model=stage_cfg.sm, solvent=self.config.general.solvent)

                main_driver = self._get_driver(
                    stage_cfg.prog,
                    stage_dir,
                    method=stage_cfg.func,
                    basis=stage_cfg.basis,
                    **opts
                )

                # 2. Create Filter via Factory
                filt = self._get_filter(stage_cfg, stage_dir, main_driver, f_type)

                # 3. Execute
                return filt.filter(ensemble, f)
        finally:
            if stage_dir.exists():
                shutil.rmtree(stage_dir)


class DefaultWorkflow(BaseWorkflow):
    def filter(self, ensemble: Ensemble) -> Ensemble:
        print(f'* Starting workflow with {len(ensemble)} conformers')

        ensemble = self._run_stage(
            '1_prescreening',
            self.config.prescreening, ensemble, 'energy' if self.config.general.gas_phase else 'gsolv')

        ensemble = self._run_stage(
            '2_screening',
            self.config.screening, ensemble, 'gibbs' if self.config.screening.gsolv_included else 'gibbs_xtb')

        return ensemble


def main():
    parser = argparse.ArgumentParser(description='pymerk Workflow Runner')
    parser.add_argument('conformers', type=pathlib.Path, help='XYZ file with conformers')
    parser.add_argument('-i', '--input', type=pathlib.Path, help='TOML config file')
    parser.add_argument('-o', '--output', type=pathlib.Path, default='final_ensemble.xyz')
    parser.add_argument('-c', '--charge', type=int, default=0)
    parser.add_argument('-m', '--multiplicity', type=int, default=1)
    args = parser.parse_args()

    print(f'* This is pymerk v{pymerk.__version__}')

    # Load Config
    config_data = {}
    if args.input:
        with args.input.open('rb') as f:
            config_data = tomllib.load(f)
    config = Config(**config_data)

    # Load Ensemble
    with args.conformers.open('r') as f:
        ensemble = Ensemble.from_multi_xyz(f, args.charge, args.multiplicity)

    # Execute Workflow
    workflow = DefaultWorkflow(pathlib.Path.cwd(), config)
    final_ensemble = workflow.filter(ensemble)

    # Save Output
    with args.output.open('w') as f:
        final_ensemble.as_multi_xyz(f)


if __name__ == '__main__':
    main()
