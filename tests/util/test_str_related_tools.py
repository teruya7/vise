# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

import pytest

from vise.util.str_related_tools import str2bool, is_str_digit, is_str_int


class TestStr2Bool:
    """Tests for the str2bool function."""

    def test_true_lowercase(self):
        """Test 'true' returns True."""
        assert str2bool("true") is True

    def test_true_uppercase(self):
        """Test 'TRUE' returns True."""
        assert str2bool("TRUE") is True

    def test_true_mixed_case(self):
        """Test 'True' returns True."""
        assert str2bool("True") is True

    def test_t_lowercase(self):
        """Test 't' returns True."""
        assert str2bool("t") is True

    def test_t_uppercase(self):
        """Test 'T' returns True."""
        assert str2bool("T") is True

    def test_false_lowercase(self):
        """Test 'false' returns False."""
        assert str2bool("false") is False

    def test_false_uppercase(self):
        """Test 'FALSE' returns False."""
        assert str2bool("FALSE") is False

    def test_false_mixed_case(self):
        """Test 'False' returns False."""
        assert str2bool("False") is False

    def test_invalid_value_raises_error(self):
        """Test invalid value raises ValueError."""
        with pytest.raises(ValueError, match="invalid truth value"):
            str2bool("yes")

    def test_empty_string_raises_error(self):
        """Test empty string raises ValueError."""
        with pytest.raises(ValueError, match="invalid truth value"):
            str2bool("")


    def test_number_string_raises_error(self):
        """Test number string raises ValueError."""
        with pytest.raises(ValueError, match="invalid truth value"):
            str2bool("1")


class TestIsStrDigit:
    """Tests for the is_str_digit function."""

    def test_integer_string(self):
        """Test integer string returns True."""
        assert is_str_digit("42") is True

    def test_negative_integer(self):
        """Test negative integer string returns True."""
        assert is_str_digit("-42") is True

    def test_float_string(self):
        """Test float string returns True."""
        assert is_str_digit("3.14") is True

    def test_negative_float(self):
        """Test negative float string returns True."""
        assert is_str_digit("-3.14") is True

    def test_scientific_notation(self):
        """Test scientific notation returns True."""
        assert is_str_digit("1e-5") is True
        assert is_str_digit("2.5E+10") is True

    def test_zero(self):
        """Test zero returns True."""
        assert is_str_digit("0") is True
        assert is_str_digit("0.0") is True

    def test_non_digit_string(self):
        """Test non-digit string returns False."""
        assert is_str_digit("hello") is False

    def test_mixed_string(self):
        """Test mixed alphanumeric string returns False."""
        assert is_str_digit("12abc") is False

    def test_empty_string(self):
        """Test empty string returns False."""
        assert is_str_digit("") is False


class TestIsStrInt:
    """Tests for the is_str_int function."""

    def test_integer_string(self):
        """Test integer string returns True."""
        assert is_str_int("42") is True

    def test_negative_integer(self):
        """Test negative integer string returns True."""
        assert is_str_int("-42") is True

    def test_zero(self):
        """Test zero returns True."""
        assert is_str_int("0") is True

    def test_float_string(self):
        """Test float string returns False."""
        assert is_str_int("3.14") is False

    def test_float_with_zero_decimal(self):
        """Test float with .0 returns True (is effectively an integer)."""
        # Note: "5.0" parses as float 5.0, but int("5.0") raises ValueError
        assert is_str_int("5.0") is False

    def test_non_digit_string(self):
        """Test non-digit string returns False."""
        assert is_str_int("hello") is False

    def test_empty_string(self):
        """Test empty string returns False."""
        assert is_str_int("") is False

    def test_large_integer(self):
        """Test large integer string returns True."""
        assert is_str_int("123456789012345678901234567890") is True
