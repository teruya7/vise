# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Mixin classes for file serialization.

This module provides abstract base classes (mixins) that add file
serialization capabilities (JSON, CSV, YAML) to dataclasses and
other data structures.
"""

import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional, TypeVar

import pandas as pd
import yaml

from vise.util.logger import get_logger

logger = get_logger(__name__)

T = TypeVar("T", bound="ToFileMixIn")


class ToFileMixIn(ABC):
    """Base mixin providing automatic filename generation.

    Generates default filenames based on class name using snake_case
    convention (e.g., MyDataClass -> my_data_class).
    """

    @property
    def _filename(self) -> str:
        """Generate snake_case filename from class name.

        Converts CamelCase class names to snake_case.
        Example: 'ClassForThis' -> 'class_for_this'

        See: https://stackoverflow.com/questions/7322028/
        """
        class_name = self.__class__.__name__
        return re.sub("(?<!^)(?=[A-Z])", "_", class_name).lower()


class ToJsonFileMixIn(ToFileMixIn, ABC):
    """Mixin for JSON file serialization."""

    def to_json_file(
        self,
        filename: Optional[str] = None,
        suffix: Optional[str] = None,
    ) -> None:
        """Write the object to a JSON file.

        Args:
            filename: Output filename. If None, uses auto-generated name.
            suffix: Optional suffix to append to the filename.
        """
        filename = filename or self.fname_w_suffix(suffix) or self.json_filename
        Path(filename).write_text(self.to_json())

    @abstractmethod
    def to_json(self) -> str:
        """Serialize the object to a JSON string.

        Returns:
            JSON-formatted string representation.
        """
        pass

    @property
    def json_filename(self) -> str:
        """Default JSON filename for this object."""
        return self._filename + ".json"

    def fname_w_suffix(self, suffix: Optional[str]) -> Optional[str]:
        """Generate filename with optional suffix.

        Args:
            suffix: Suffix to append before extension.

        Returns:
            Filename with suffix, or None if no suffix provided.
        """
        if suffix:
            return f"{self._filename}_{suffix}.json"
        return None


class ToCsvFileMixIn(ToFileMixIn, ABC):
    """Mixin for CSV file serialization using pandas DataFrame."""

    def to_csv_file(self, filename: Optional[str] = None) -> None:
        """Write the object to a CSV file.

        Args:
            filename: Output filename. If None, uses auto-generated name.
        """
        self.to_dataframe.to_csv(filename or self._csv_filename, index=False)

    @classmethod
    def from_csv_file(cls: type[T], filename: str) -> T:
        """Load an object from a CSV file.

        Args:
            filename: Path to the CSV file.

        Returns:
            New instance of the class.
        """
        return cls.from_dataframe(pd.read_csv(filename))

    @property
    @abstractmethod
    def to_dataframe(self) -> pd.DataFrame:
        """Convert the object to a pandas DataFrame.

        Returns:
            DataFrame representation of the object.
        """
        pass

    @classmethod
    @abstractmethod
    def from_dataframe(cls: type[T], df: pd.DataFrame) -> T:
        """Create an instance from a pandas DataFrame.

        Args:
            df: Source DataFrame.

        Returns:
            New instance of the class.
        """
        pass

    @property
    def _csv_filename(self) -> str:
        """Default CSV filename for this object."""
        return self._filename + ".csv"


class ToYamlFileMixIn(ToFileMixIn, ABC):
    """Mixin for YAML file serialization."""

    def to_yaml_file(self, filename: Optional[str] = None) -> None:
        """Write the object to a YAML file.

        Args:
            filename: Output filename. If None, uses auto-generated name.
        """
        filename = filename or self._yaml_filename
        Path(filename).write_text(self.to_yaml())

    def to_yaml(self) -> str:
        """Serialize the object to a YAML string.

        Returns:
            YAML-formatted string representation.
        """
        return yaml.dump(self.as_dict())

    def as_dict(self) -> Optional[Dict[str, Any]]:
        """Convert the object to a dictionary.

        Override this method to provide custom serialization.
        Not abstract because subclasses may override to_yaml/from_yaml instead.

        Returns:
            Dictionary representation, or None if not implemented.
        """
        return None

    @classmethod
    def from_yaml(cls: type[T], filename: str) -> T:
        """Load an object from a YAML file.

        Args:
            filename: Path to the YAML file.

        Returns:
            New instance of the class.
        """
        with open(filename) as file:
            data = yaml.safe_load(file)
        if hasattr(cls, "from_dict"):
            return cls.from_dict(data)
        return cls(**data)

    @property
    def _yaml_filename(self) -> str:
        """Default YAML filename for this object."""
        return self._filename + ".yaml"
