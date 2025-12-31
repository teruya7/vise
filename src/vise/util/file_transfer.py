# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""File transfer utilities for vise.

This module provides classes for transferring files between directories
using different methods (move, copy, symbolic link) with support for
compressed files.
"""

import os
import re
import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Type

from monty.io import zopen

from vise.error import ViseError
from vise.util.logger import get_logger

logger = get_logger(__name__)


class FileTransfer(ABC):
    """Abstract base class for file transfer operations.

    Provides a common interface for different file transfer methods.

    Attributes:
        abs_in_file: Absolute path to the source file.
    """

    def __init__(self, abs_file_path: Path) -> None:
        """Initialize with the source file path.

        Args:
            abs_file_path: Absolute path to the file to transfer.
        """
        self.abs_in_file = abs_file_path

    @property
    def file_name(self) -> str:
        """Get the file name from the path."""
        return Path(self.abs_in_file).name

    @abstractmethod
    def transfer(self, abs_output_dir: Path) -> None:
        """Transfer the file to the output directory.

        Args:
            abs_output_dir: Absolute path to the destination directory.
        """
        pass

    def to(self, abs_output_dir: Path) -> Path:
        """Get the destination path for this file.

        Args:
            abs_output_dir: Absolute path to the destination directory.

        Returns:
            Full path where the file will be placed.
        """
        return abs_output_dir / self.file_name


class FileMove(FileTransfer):
    """Transfer file by moving it."""

    def transfer(self, abs_output_dir: Path) -> None:
        """Move the file to the output directory."""
        shutil.move(self.abs_in_file, self.to(abs_output_dir))


class FileCopy(FileTransfer):
    """Transfer file by copying it (supports compressed files)."""

    def transfer(self, abs_output_dir: Path) -> None:
        """Copy the file to the output directory.

        Uses monty.io.zopen to support compressed files.
        """
        with zopen(self.abs_in_file, "rb") as fin, \
                zopen(self.to(abs_output_dir), "wb") as fout:
            shutil.copyfileobj(fin, fout)


class FileLink(FileTransfer):
    """Transfer file by creating a symbolic link."""

    def transfer(self, abs_output_dir: Path) -> None:
        """Create a symbolic link to the file in the output directory."""
        os.symlink(self.abs_in_file, self.to(abs_output_dir))


def transfer_instance(
    transfer_type: str, filename: Path
) -> FileTransfer:
    """Create a FileTransfer instance based on the transfer type.

    Args:
        transfer_type: Type of transfer - 'move', 'copy', or 'link'
                      (only first character is checked).
        filename: Path to the file to transfer.

    Returns:
        Appropriate FileTransfer subclass instance.

    Raises:
        ViseFileTransferError: If transfer type is not recognized.
    """
    transfer_classes: Dict[str, Type[FileTransfer]] = {
        "m": FileMove,
        "c": FileCopy,
        "l": FileLink,
    }

    initial = transfer_type[0].lower()
    if initial not in transfer_classes:
        raise ViseFileTransferError(
            f"Transfer type '{initial}' is not supported. "
            f"Choose from: m (move), c (copy), or l (link)."
        )

    return transfer_classes[initial](filename)


class FileTransfers:
    """Collection of file transfers.

    Manages multiple file transfer operations based on a configuration
    dictionary specifying transfer types for each file.
    """

    def __init__(
        self, file_transfer_types: Dict[str, str], path: Path
    ) -> None:
        """Initialize with file transfer specifications.

        Args:
            file_transfer_types: Dict mapping filename to transfer type.
            path: Base path where source files are located.
        """
        file_transfers: List[FileTransfer] = []

        for filename, transfer_type in file_transfer_types.items():
            file_path = path.absolute() / filename

            if not file_path.is_file():
                logger.warning(f"{file_path} does not exist.")
            elif file_path.stat().st_size == 0:
                logger.warning(f"{file_path} is empty.")
            else:
                file_transfers.append(
                    transfer_instance(transfer_type, file_path)
                )

        self.file_transfers = file_transfers

    def delete_file_transfers(self, keywords: List[str]) -> None:
        """Remove files matching any of the keywords from transfer list.

        Args:
            keywords: List of patterns to match against filenames.
        """
        pattern = re.compile("|".join(keywords))
        for file_transfer in list(self.file_transfers):
            if pattern.search(file_transfer.file_name):
                self.file_transfers.remove(file_transfer)

    def transfer(self, output_dir: Path = Path(".")) -> None:
        """Execute all file transfers to the output directory.

        Args:
            output_dir: Destination directory (default: current directory).
        """
        for transfer in self.file_transfers:
            transfer.transfer(output_dir.absolute())


class ViseFileTransferError(ViseError):
    """Exception for file transfer errors."""

    pass
