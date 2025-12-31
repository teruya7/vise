# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""KPOINTS file generation for VASP.

This module provides ViseKpoints, an enhanced version of pymatgen's Kpoints
class with improved string formatting for k-point output.
"""

from typing import List

from pymatgen.io.vasp.sets import Kpoints


class ViseKpoints(Kpoints):
    """Enhanced Kpoints class with improved string formatting.

    Extends pymatgen's Kpoints to provide precise formatting of k-point
    coordinates, particularly for line-mode band structure calculations.
    """

    def __str__(self) -> str:
        """Format the KPOINTS file content.

        Returns:
            Properly formatted KPOINTS file content as a string.

        Note:
            This method improves on pymatgen's formatting by:
            - Using high precision for multi-point k-paths
            - Properly formatting weights and labels
            - Handling line-mode with path segment markers
        """
        lines: List[str] = [
            self.comment,
            str(self.num_kpts),
            self.style.name,
        ]

        style_char = self.style.name.lower()[0]
        is_line_mode = style_char == "l"

        if is_line_mode:
            lines.append(self.coord_type)

        for i, kpt in enumerate(self.kpts):
            # Format k-point coordinates
            if len(self.kpts) == 1:
                # Single point: simple formatting
                kpt_str = " ".join(str(x) for x in kpt)
            else:
                # Multiple points: high precision
                kpt_str = " ".join(f"{x:20.17f}" for x in kpt)

            if is_line_mode:
                # Line mode: add label
                kpt_str += f" ! {self.labels[i]}"
                # Add blank line between path segments
                if i % 2 == 1:
                    kpt_str += "\n"
            elif self.num_kpts > 0:
                # Explicit k-points: add weight and optional label
                if self.labels is not None:
                    kpt_str += f"{self.kpts_weights[i]:4} {self.labels[i]}"
                else:
                    kpt_str += f"{self.kpts_weights[i]:4}"

            lines.append(kpt_str)

        # Add shift for automatic k-point meshes if non-zero
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            shift_str = " ".join(str(x) for x in self.kpts_shift)
            lines.append(shift_str)

        return "\n".join(lines) + "\n"
