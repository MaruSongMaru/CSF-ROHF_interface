#!/usr/bin/env python3
"""CSF-specific helper functions for MOLCAS validation."""

from __future__ import annotations

import re
from typing import Iterable


def calculate_csf_active_electrons(csf_stepvec: Iterable[int]) -> int:
    """Calculate the active electron count implied by a CSF step vector."""
    electron_map = {
        0: 0,
        1: 1,
        2: 1,
        3: 2,
    }

    active_electrons = 0
    for step in csf_stepvec:
        try:
            step_value = int(step)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid CSF step vector element: {step!r}") from exc

        if step_value not in electron_map:
            raise ValueError(f"Unsupported CSF step vector element: {step_value}")

        active_electrons += electron_map[step_value]

    return active_electrons


def calculate_csf_spin_twice(csf_stepvec: Iterable[int]) -> int:
    """Calculate 2S implied by a CSF step vector.

    The cumulative value must never become negative while scanning left to right.
    """
    spin_map = {
        0: 0,
        1: 1,
        2: -1,
        3: 0,
    }

    spin_twice = 0
    for step in csf_stepvec:
        try:
            step_value = int(step)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid CSF step vector element: {step!r}") from exc

        if step_value not in spin_map:
            raise ValueError(f"Unsupported CSF step vector element: {step_value}")

        spin_twice += spin_map[step_value]
        if spin_twice < 0:
            raise ValueError(
                f"Invalid CSF step vector: cumulative 2S became negative ({spin_twice})"
            )

    return spin_twice


def extract_active_orbitals(log_content: str) -> int:
    """Extract the active orbital count from MOLCAS log content."""
    matches = re.findall(r"Number of active orbitals\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'Number of active orbitals' in MOLCAS log")

    return int(matches[-1])


def extract_active_shell_electrons(log_content: str) -> int:
    """Extract the active-shell electron count from MOLCAS log content."""
    matches = re.findall(r"Number of electrons in active shells\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'Number of electrons in active shells' in MOLCAS log")

    return int(matches[-1])


def extract_state_symmetry(log_content: str) -> int:
    """Extract the state symmetry / multiplicity value from MOLCAS log content."""
    matches = re.findall(r"State symmetry\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'State symmetry' in MOLCAS log")

    return int(matches[-1])


def validate_csf_against_molcas(log_content: str, csf_stepvec: Iterable[int]) -> None:
    """Validate the CSF against the first-RASSCF MOLCAS log information."""
    csf_stepvec_list = list(csf_stepvec)

    active_orbitals = extract_active_orbitals(log_content)
    if len(csf_stepvec_list) != active_orbitals:
        raise ValueError(
            f"CSF step vector length ({len(csf_stepvec_list)}) does not match the number of active orbitals "
            f"({active_orbitals}) at first RASSCF iteration"
        )

    active_shell_electrons = extract_active_shell_electrons(log_content)
    csf_active_electrons = calculate_csf_active_electrons(csf_stepvec_list)
    if csf_active_electrons != active_shell_electrons:
        raise ValueError(
            f"CSF active-electron count ({csf_active_electrons}) does not match the number of electrons in active shells "
            f"({active_shell_electrons}) at first RASSCF iteration"
        )

    state_symmetry = extract_state_symmetry(log_content)
    csf_spin_twice = calculate_csf_spin_twice(csf_stepvec_list)
    csf_multiplicity = csf_spin_twice + 1
    if csf_multiplicity != state_symmetry:
        raise ValueError(
            f"CSF multiplicity ({csf_multiplicity}) does not match MOLCAS state symmetry ({state_symmetry}) at first RASSCF iteration"
        )