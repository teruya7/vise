# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Matplotlib utilities for vise.

This module provides custom formatters and utilities for matplotlib plots,
particularly for formatting axis tick labels.
"""

from typing import Union

from matplotlib.ticker import FuncFormatter


def format_tick_value(tick_value: float, pos: int) -> Union[int, float, str]:
    """Format tick values, converting integer-valued floats to integers.

    This formatter converts values like 0.0 to 0 in plots for cleaner display.
    The `pos` parameter is required by matplotlib's FuncFormatter interface
    but is not used in this implementation.

    Args:
        tick_value: The tick value to format.
        pos: Position of the tick (unused, required by FuncFormatter).

    Returns:
        Formatted tick value: int if the float is effectively an integer,
        otherwise the rounded float or string representation.

    Examples:
        >>> format_tick_value(0.0, 0)
        0
        >>> format_tick_value(2.5, 0)
        2.5

    Note:
        To use this formatter with a matplotlib axis:
            ax.xaxis.set_major_formatter(float_to_int_formatter)
            ax.yaxis.set_major_formatter(float_to_int_formatter)
    """
    if isinstance(tick_value, float):
        rounded_value = round(tick_value, ndigits=10)
        if rounded_value.is_integer():
            return int(rounded_value)
        else:
            return rounded_value
    else:
        return str(tick_value)


# Pre-configured formatter instance for convenience
float_to_int_formatter = FuncFormatter(format_tick_value)

# Backward compatibility alias
my_formatter = format_tick_value
