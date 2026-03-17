"""Tests for pure functions in pdbe_sifts.base.utils."""

import os
from unittest.mock import patch

from pdbe_sifts.base.utils import SiftsAction, get_allocated_cpus, get_cpu_count


class TestGetAllocatedCpus:
    def test_slurm_env_var(self):
        with patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "8"}):
            assert get_allocated_cpus() == 8

    def test_fallback_returns_positive_int(self):
        env = {
            k: v for k, v in os.environ.items() if k != "SLURM_CPUS_PER_TASK"
        }
        with patch.dict(os.environ, env, clear=True):
            result = get_allocated_cpus()
            assert isinstance(result, int)
            assert result >= 1


class TestGetCpuCount:
    def test_sifts_n_proc_env_var(self):
        with patch.dict(os.environ, {"SIFTS_N_PROC": "4"}):
            assert get_cpu_count() == 4

    def test_fallback_to_allocated_cpus(self):
        env = {k: v for k, v in os.environ.items() if k != "SIFTS_N_PROC"}
        with patch.dict(os.environ, env, clear=True):
            result = get_cpu_count()
            assert isinstance(result, int)
            assert result >= 1

    def test_sifts_n_proc_overrides_slurm(self):
        with patch.dict(
            os.environ,
            {"SIFTS_N_PROC": "2", "SLURM_CPUS_PER_TASK": "16"},
        ):
            assert get_cpu_count() == 2


class TestSiftsAction:
    """SiftsAction extends argparse.Action — test default resolution logic."""

    def _make_action(self, envvar=None, confvar=None, default=None):
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--foo",
            action=SiftsAction,
            envvar=envvar,
            confvar=confvar,
            default=default,
        )
        return parser

    def test_explicit_default(self):
        p = self._make_action(default="explicit")
        args = p.parse_args([])
        assert args.foo == "explicit"

    def test_env_var_sets_default(self):
        with patch.dict(os.environ, {"MY_VAR": "from_env"}):
            p = self._make_action(envvar="MY_VAR")
            args = p.parse_args([])
            assert args.foo == "from_env"

    def test_confvar_fallback(self):
        env = {k: v for k, v in os.environ.items() if k != "MY_VAR"}
        with patch.dict(os.environ, env, clear=True):
            p = self._make_action(envvar="MY_VAR", confvar="from_conf")
            args = p.parse_args([])
            assert args.foo == "from_conf"

    def test_cli_value_overrides_default(self):
        p = self._make_action(default="default_val")
        args = p.parse_args(["--foo", "cli_val"])
        assert args.foo == "cli_val"
