import argparse
import pathlib
import shutil
from typing import Any, Callable

import pymerk
from pymerk.driver import XtbDriver, VlxDriver, BaseDriver, OrcaDriver
from pymerk.ensemble import Ensemble
from pymerk.scripts import Config
from pymerk.scripts.filter import EnergyFilter, SelectDriver, OptFilter, MacroOptFilter, BoltzmannFilter

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
            opts = {
                'imagthr': self.config.general.imagthr,
                'sthr': self.config.general.sthr,
            }
            if not self.config.general.gas_phase:
                opts.update({
                    'solvatation_model': self.config.general.sm_rrho,
                    'solvent': self.config.general.solvent
                })
            opts.update(overrides)
            return XtbDriver(stage_dir, paths.xtb, **opts)

        if prog == 'vlx':
            if not paths.vlx:
                raise RuntimeError('VeloxChem path not configured.')
            return VlxDriver(stage_dir, paths.vlx, **overrides)

        if prog == 'orca':
            if not paths.orca:
                raise RuntimeError('VeloxChem path not configured.')
            return OrcaDriver(stage_dir, paths.orca, nprocs=self.config.paths.orca_nprocs, **overrides)

        raise ValueError(f'Unsupported driver: {prog}')

    def _execute_stage(self, name: str, stage_func: Callable[[pathlib.Path, Any], Ensemble]) -> Ensemble:
        """
        Generic wrapper for any stage type
        Handles directory creation, header printing, and cleanup.
        """
        self._print_header(f'Stage: {name}')
        stage_dir = self.workdir / name
        log_file = self.workdir / f'{name}.log'
        xyz_file = self.workdir / f'{name}.selected.xyz'
        stage_dir.mkdir(exist_ok=True)

        try:
            with log_file.open('w') as f:
                ensemble = stage_func(stage_dir, f)
            with xyz_file.open('w') as f:
                ensemble.as_multi_xyz(f)
            return ensemble
        finally:
            if stage_dir.exists():
                shutil.rmtree(stage_dir)


class DefaultWorkflow(BaseWorkflow):

    def _resolve_filter_components(self, stage_name: str, stage_cfg: Any):
        """Helper to resolve component logic based on general gas-phase settings."""
        if self.config.general.gas_phase:
            return SelectDriver.NONE, SelectDriver.NONE, 'ΔE'

        if stage_name in ['2_screening', '4_refinement']:
            gsolv = SelectDriver.MAIN if getattr(stage_cfg, 'gsolv_included', False) else (
                SelectDriver.AUX if not self.config.general.gas_phase else SelectDriver.NONE
            )
            gtrv = SelectDriver.AUX if self.config.general.evaluate_rrho else SelectDriver.MAIN
            label = 'G' if gsolv == SelectDriver.NONE else 'G*'
        else:
            # Prescreening defaults
            gsolv = SelectDriver.AUX if not self.config.general.gas_phase else SelectDriver.NONE
            gtrv = SelectDriver.NONE
            label = 'E' if gsolv == SelectDriver.NONE else 'g*'

        return gsolv, gtrv, label

    def _run_energy_filter(
            self,
            stage_dir: pathlib.Path,
            log_output: Any,
            ensemble: Ensemble,
            stage_cfg: Any,
            gsolv: SelectDriver,
            gtrv: SelectDriver,
            label: str
    ) -> Ensemble:
        """Specific logic for setting up an Energy Filter."""

        # 1. Setup main driver options
        opts = {'method': stage_cfg.func, 'basis': getattr(stage_cfg, 'basis', None)}

        # Determine if the main driver should handle solvation
        if hasattr(stage_cfg, 'sm') and stage_cfg.sm != '' and getattr(stage_cfg, 'gsolv_included', True):
            solvent = self.config.general.solvent
            if getattr(stage_cfg, 'alternate_solvent', None) is not None:
                solvent = stage_cfg.alternate_solvent
            opts.update(solvatation_model=stage_cfg.sm, solvent=solvent)

        main_driver = self._get_driver(stage_cfg.prog, stage_dir, **opts)

        # 2. Setup auxiliary driver (if needed)
        aux_driver = None
        if gsolv == SelectDriver.AUX or gtrv == SelectDriver.AUX:
            aux_driver = self._get_driver('xtb', stage_dir, version=stage_cfg.gfnv)

        # 3. Apply filter
        filt = EnergyFilter(
            main_driver, stage_cfg.threshold / AU_TO_KCAL,
            gsolv, gtrv, aux_driver=aux_driver, label=label
        )
        return filt.filter(ensemble, log_output)

    def _run_boltzmann_filter(
            self,
            stage_dir: pathlib.Path,
            log_output: Any,
            ensemble: Ensemble,
            stage_cfg: Any,
            gsolv: SelectDriver,
            gtrv: SelectDriver,
            label: str
    ) -> Ensemble:
        """Specific logic for setting up an Energy Filter."""

        # 1. Setup main driver options
        opts = {'method': stage_cfg.func, 'basis': getattr(stage_cfg, 'basis', None)}

        # Determine if the main driver should handle solvation
        if hasattr(stage_cfg, 'sm') and stage_cfg.sm != '' and getattr(stage_cfg, 'gsolv_included', True):
            solvent = self.config.general.solvent
            if getattr(stage_cfg, 'alternate_solvent', None) is not None:
                solvent = stage_cfg.alternate_solvent
            opts.update(solvatation_model=stage_cfg.sm, solvent=solvent)

        main_driver = self._get_driver(stage_cfg.prog, stage_dir, **opts)

        # 2. Setup auxiliary driver (if needed)
        aux_driver = None
        if gsolv == SelectDriver.AUX or gtrv == SelectDriver.AUX:
            aux_driver = self._get_driver('xtb', stage_dir, version=stage_cfg.gfnv)

        # 3. Apply filter
        filt = BoltzmannFilter(
            main_driver, stage_cfg.threshold,
            gsolv, gtrv, aux_driver=aux_driver, label=label
        )
        return filt.filter(ensemble, log_output)

    def _run_opt_filter(
            self,
            stage_dir: pathlib.Path,
            log_output: Any,
            ensemble: Ensemble,
            stage_cfg: Any,
            use_solvent: bool,
            gtrv: SelectDriver,
            label: str
    ) -> Ensemble:
        """Specific logic for setting up an Opt Filter."""

        # 1. Setup main driver options
        opts = {'method': stage_cfg.func, 'basis': getattr(stage_cfg, 'basis', None)}

        # Determine if the main driver should handle solvation
        if hasattr(stage_cfg, 'sm') and stage_cfg.sm != '' and getattr(stage_cfg, 'gsolv_included', True):
            solvent = self.config.general.solvent
            if getattr(stage_cfg, 'alternate_solvent', None) is not None:
                solvent = stage_cfg.alternate_solvent
            opts.update(solvatation_model=stage_cfg.sm, solvent=solvent)

        main_driver = self._get_driver(stage_cfg.prog, stage_dir, **opts)

        # 2. Setup auxiliary driver (if needed)
        aux_driver = None
        if gtrv == SelectDriver.AUX:
            aux_driver = self._get_driver('xtb', stage_dir, version=stage_cfg.gfnv)

        # 3. Apply filter
        if getattr(stage_cfg, 'macrocycles', False):
            filt = MacroOptFilter(
                main_driver, stage_cfg.threshold / AU_TO_KCAL,
                use_solvent, gtrv, aux_driver=aux_driver, label=label,
                maxcycles=stage_cfg.maxcyc, optcycles=stage_cfg.optcycles, gradthr=stage_cfg.gradthr,
                optlevel=stage_cfg.optlevel
            )
        else:
            filt = OptFilter(
                main_driver, stage_cfg.threshold / AU_TO_KCAL,
                use_solvent, gtrv, aux_driver=aux_driver, label=label,
                maxcycles=stage_cfg.maxcyc, optlevel=stage_cfg.optlevel
            )
        return filt.filter(ensemble, log_output)

    def filter(self, ensemble: Ensemble) -> Ensemble:
        print(f'* Starting workflow with {len(ensemble)} conformers')
        print(f'* Workdir: {self.workdir}')

        # 1. Prescreening
        g, t, label = self._resolve_filter_components('1_prescreening', self.config.prescreening)
        ensemble = self._execute_stage(
            '1_prescreening',
            lambda d, f: self._run_energy_filter(d, f, ensemble, self.config.prescreening, g, t, label)
        )

        # 2. Screening
        g, t, label = self._resolve_filter_components('2_screening', self.config.screening)
        ensemble = self._execute_stage(
            '2_screening',
            lambda d, f: self._run_energy_filter(d, f, ensemble, self.config.screening, g, t, label)
        )

        # 3. Optimize
        t = SelectDriver.AUX if self.config.general.evaluate_rrho else SelectDriver.MAIN
        label = 'G' if self.config.general.gas_phase else 'G*'

        ensemble = self._execute_stage(
            '3_optimize',
            lambda d, f: self._run_opt_filter(
                d, f, ensemble, self.config.optimization, not self.config.general.gas_phase, t, label)
        )

        # 4. Refinement
        g, t, label = self._resolve_filter_components('4_refinement', self.config.refinement)
        ensemble = self._execute_stage(
            '4_refinement',
            lambda d, f: self._run_boltzmann_filter(d, f, ensemble, self.config.refinement, g, t, label)
        )

        return ensemble


def main():
    parser = argparse.ArgumentParser(description='pymerk Workflow Runner')
    parser.add_argument('conformers', type=pathlib.Path, help='XYZ file with conformers')
    parser.add_argument('-i', '--input', type=pathlib.Path, help='TOML config file')
    parser.add_argument('-o', '--output', type=pathlib.Path, default='final_ensemble.xyz')
    parser.add_argument('-c', '--charge', type=int, default=0)
    parser.add_argument('-m', '--multiplicity', type=int, default=1)
    parser.add_argument('-w', '--workdir', type=pathlib.Path, default=pathlib.Path.cwd(), help='Working directory')

    args = parser.parse_args()

    print(f'* This is pymerk v{pymerk.__version__}')

    # Load Config
    if args.input:
        with args.input.open('rb') as f:
            config = Config.from_toml(f)
    else:
        config = Config()

    # Load Ensemble
    with args.conformers.open('r') as f:
        ensemble = Ensemble.from_multi_xyz(f, args.charge, args.multiplicity)

    # Execute Workflow
    workflow = DefaultWorkflow(args.workdir, config)
    final_ensemble = workflow.filter(ensemble)

    # Save Output
    with args.output.open('w') as f:
        final_ensemble.as_multi_xyz(f)


if __name__ == '__main__':
    main()
