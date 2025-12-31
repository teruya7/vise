# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Extended Enum base class for vise.

This module provides an enhanced Enum class with additional utility methods
for string conversion and member listing.
"""

from enum import Enum, unique
from typing import List, TypeVar

T = TypeVar("T", bound="ExtendedEnum")


@unique
class ExtendedEnum(Enum):
    """Enhanced Enum with additional utility methods.

    Provides string conversion, member lookup by name, and
    convenient methods for listing enum members and values.
    """

    def __repr__(self) -> str:
        """Return the enum member name."""
        return self.name

    def __str__(self) -> str:
        """Return the enum member name as a string."""
        return self.name

    @classmethod
    def from_string(cls: type[T], name: str) -> T:
        """Create an enum member from its string name.

        Args:
            name: The name of the enum member.

        Returns:
            The matching enum member.

        Raises:
            AttributeError: If no member with the given name exists.

        Examples:
            >>> class Color(ExtendedEnum):
            ...     RED = 1
            >>> Color.from_string("RED")
            RED
        """
        for member in cls:
            if member.name == name:
                return member
        raise AttributeError(
            f"{name!r} is not a valid member of {cls.__name__}. "
            f"Valid names: {cls.names_string()}"
        )

    @classmethod
    def names_string(cls) -> str:
        """Return a comma-separated string of all member names.

        Returns:
            String containing all member names separated by commas.
        """
        return ", ".join(member.name for member in cls)

    @classmethod
    def name_list(cls: type[T]) -> List[T]:
        """Return a list of all enum members.

        Returns:
            List containing all enum members.
        """
        return list(cls)

    @classmethod
    def name_str_list(cls) -> List[str]:
        """Return a list of all member names as strings.

        Returns:
            List of member name strings.
        """
        return [member.name for member in cls]

    @classmethod
    def values(cls) -> list:
        """Return a list of all member values.

        Returns:
            List containing the value of each enum member.
        """
        return [member.value for member in cls]
