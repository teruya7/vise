# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""User settings loader for VISE configuration.

This module provides the UserSettings class that discovers and loads
user configuration files (YAML) from the directory hierarchy. Settings
files are searched from the current working directory up to the root,
with files closer to the current directory taking precedence.

Configuration files can be named either:
- vise.yaml (visible)
- .vise.yaml (hidden)
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from vise.util.logger import get_logger

logger = get_logger(__name__)


class UserSettings:
    """Loader for user configuration files from the directory hierarchy.

    This class searches for YAML configuration files starting from the
    current working directory and walking up to the root directory.
    Files are loaded in order from root to current directory, so that
    settings in directories closer to the current directory override
    those in parent directories.

    Attributes:
        yaml_filename: The base filename to search for.
        yaml_files_from_root_dir: List of found YAML files, ordered from
            root to current directory.
    """

    def __init__(self, yaml_filename: str) -> None:
        """Initialize and search for YAML configuration files.

        Args:
            yaml_filename: Base filename to search for (e.g., 'vise.yaml').
                Also searches for hidden variant (e.g., '.vise.yaml').
        """
        self._cwd = Path.cwd()
        self.yaml_filename = yaml_filename
        self.yaml_files_from_root_dir = self._make_yaml_file_list()

    def _make_yaml_file_list(self) -> List[Path]:
        """Search for YAML files from current directory to root.

        Returns:
            List of found YAML file paths, ordered from root to current
            directory (files found later in the list override earlier ones).
        """
        result: List[Path] = []

        current_dir = self._cwd
        while True:
            # Search for both visible and hidden variants
            filenames = [self.yaml_filename, "." + self.yaml_filename]
            file_paths = [current_dir / filename for filename in filenames]

            for file_path in file_paths:
                if file_path.exists():
                    result.append(file_path)

            if current_dir == Path("/"):
                break
            else:
                current_dir = current_dir.parent

        # Reverse to get root-to-current order
        return list(reversed(result))

    @property
    def user_settings(self) -> Dict[str, Any]:
        """Load and merge settings from all found YAML files.

        Files are processed in order from root to current directory,
        so settings in files closer to the current directory override
        those in parent directories. Dictionary values are merged rather
        than replaced.

        Path-like string values (containing '/') are resolved relative
        to the containing YAML file's directory.

        Returns:
            Merged dictionary of all user settings.
        """
        result: Dict[str, Any] = {}

        for file_path in self.yaml_files_from_root_dir:
            logger.info(f"Setting file: {file_path} is parsed...")
            try:
                with open(str(file_path), "r") as fin:
                    settings = yaml.load(fin, Loader=yaml.SafeLoader)
                    for key, value in settings.items():
                        if key in result:
                            logger.info(
                                f"key {key} was overridden to {value} by {file_path}"
                            )
                            # Merge dictionaries instead of replacing
                            if isinstance(value, dict):
                                result[key].update(value)
                                continue

                        # Resolve path-like values relative to YAML file
                        if self.is_path(value):
                            value = file_path.parent / value
                        result[key] = value
            except AttributeError:
                # Skip files with invalid/empty content
                pass

        logger.info(f"-- Settings from {self.yaml_filename}:")
        logger.info(", ".join(f"{k}: {v} " for k, v in result.items()))
        return result

    @staticmethod
    def is_path(value: Any) -> bool:
        """Check if a value appears to be a file path.

        Args:
            value: Value to check.

        Returns:
            True if value is a string containing at least one '/'.
        """
        return isinstance(value, str) and bool(re.match(r"\S*/\S*", value))
