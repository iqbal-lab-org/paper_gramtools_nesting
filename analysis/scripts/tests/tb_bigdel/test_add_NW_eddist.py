from unittest.mock import patch

import pytest
from pysam import AlignedSegment

from tb_bigdel.add_NW_eddist import (
    is_primary_alignment,
    get_ref_pos,
    count_matches,
    get_NW_edit_distance,
)


@patch("pysam.AlignedSegment", autospec=True, spec_set=True)
class TestSamFlags:
    def test_bit_flag_0_is_primary_alignment(self, mock_read):
        mock_read.flag = 0
        assert is_primary_alignment(mock_read)

    def test_bit_flag_256_is_not_primary_alignment(self, mock_read):
        mock_read.flag = 256
        assert not is_primary_alignment(mock_read)


@patch("pysam.AlignedSegment")
class TestRefPos:
    def test_first_is_ref_consuming(self, mock_read):
        mock_read.cigartuples = [(0, 20), (4, 20)]
        mock_read.pos = 20
        assert get_ref_pos(mock_read) == 19

    def test_first_is_not_ref_consuming(self, mock_read):
        mock_read.cigartuples = [(4, 20), (0, 20)]
        mock_read.pos = 50
        assert get_ref_pos(mock_read) == 29


@patch("pysam.AlignedSegment")
class TestCountMatches:
    def test_no_match_return_0(self, mock_read):
        mock_read.cigartuples = [(1, 20), (2, 20)]
        assert count_matches(mock_read) == 0

    def test_matches_get_counted(self, mock_read):
        mock_read.cigartuples = [(0, 20), (7, 20), (8, 2)]
        assert count_matches(mock_read) == 42


class TestNeedlemanWunsch:
    def test_identical_sequences_returns_0(self):
        query = "ATCGT"
        target = query
        assert get_NW_edit_distance(query, target) == 0

    def test_different_sequences_returns_eddist(self):
        query = "ATCGT"
        target = "AAAAAT"
        assert get_NW_edit_distance(query, target) == 4
