from unittest.mock import patch

import pytest

from pdbe_sifts.global_mappings.scoring_function_helper import get_tax_weight

_PATCH = "pdbe_sifts.global_mappings.scoring_function_helper.NCBITaxa"


def test_same_taxid():
    assert get_tax_weight(10665, 10665) == 200


def test_child_parent(mock_ncbi):
    with patch(_PATCH, return_value=mock_ncbi):
        assert get_tax_weight(10665, 697290) == 100


def test_sibling(mock_ncbi):
    with patch(_PATCH, return_value=mock_ncbi):
        assert get_tax_weight(10665, 2750851) == 50


def test_cousin(mock_ncbi):
    with patch(_PATCH, return_value=mock_ncbi):
        assert get_tax_weight(10665, 2696339) == 25


def test_distant(mock_ncbi):
    with patch(_PATCH, return_value=mock_ncbi):
        assert get_tax_weight(10665, 3044333) == 0
