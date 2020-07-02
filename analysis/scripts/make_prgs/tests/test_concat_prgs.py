from io import StringIO, BufferedIOBase
from unittest import mock
from pathlib import Path

import pytest

from starting_prg.concat_prgs import (
    load_prg_names,
    get_file_names,
    ENDIANNESS,
    BYTES_PER_INT,
    to_bytes,
    PRGAggregator,
    get_aggregated_prgs,
    PRGDecodeError,
    PRGAggregationError,
)


class TestGetPrgs:
    @classmethod
    def setup_class(cls):
        cls.features = StringIO(
            "ref1\t1\t10\tprg_1\n" "ref1\t20\t30\tprg_2\n" "ref1\t50\t80\tprg_3\n"
        )
        cls.prg_names = ["prg_1", "prg_2", "prg_3"]

    def test_get_prg_names(self):
        result = load_prg_names(self.features)
        assert result == self.prg_names

    def test_get_file_names_incomplete_fails(self):
        with mock.patch.object(Path, "iterdir") as mocked:
            mocked.return_value = [
                Path("prg_1.bin"),
                Path("prg_3.bin"),
            ]  # Missing 'prg_2' file
            with pytest.raises(FileNotFoundError):
                result = get_file_names(Path("."), self.prg_names)

    def test_get_file_names_complete(self):
        return_paths = [
            Path("prg_1.full.bin"),
            Path("prg_2.bin"),
            Path("prg_3.ext"),
        ]

        with mock.patch.object(Path, "iterdir") as mocked:
            mocked.return_value = return_paths
            result = get_file_names(Path("dummy"), self.prg_names)

        expected = {
            "prg_1": return_paths[0],
            "prg_2": return_paths[1],
            "prg_3": return_paths[2],
        }
        assert result == expected

    def test_get_file_names_keys_ordered_by_prg_name_order(self):
        return_paths = [
            Path("prg_3.bin"),
            Path("prg_2.bin"),
            Path("prg_1.bin"),
        ]

        with mock.patch.object(Path, "iterdir") as mocked:
            mocked.return_value = return_paths
            result = get_file_names(Path("dummy"), self.prg_names)

        assert list(result.keys()) == self.prg_names


class TestByteConversions:
    def test_constants(self):
        assert BYTES_PER_INT == 4
        assert ENDIANNESS == "little"

    def test_to_bytes(self):
        int_vec = [1, 4, 5]
        result = to_bytes(int_vec)
        assert result == bytes(
            b"\x01\x00\x00\x00" b"\x04\x00\x00\x00" b"\x05\x00\x00\x00"
        )


@mock.patch.object(BufferedIOBase, "read")
@mock.patch.object(Path, "open", side_effect=lambda _: BufferedIOBase())
class TestPRGAggregation:
    def test_prg_with_invalid_integer_fails(self, op, re):
        invalid_prg = [0, 1]
        re.side_effect = [to_bytes(invalid_prg)]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        with pytest.raises(PRGDecodeError):
            result = get_aggregated_prgs(agg, file_names)

    def test_one_prg_no_variants(self, op, re):
        prg_1 = [1, 2, 1, 1, 2, 4]
        re.side_effect = [to_bytes(prg_1)]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        assert result == prg_1
        assert agg.next_allocated == 5

    def test_one_prg_with_variants(self, op, re):
        prg_1 = [1, 2, 5, 1, 6, 4, 6, 1, 7, 3, 3, 8, 4, 8]
        re.side_effect = [to_bytes(prg_1)]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        assert result == prg_1
        assert agg.next_allocated == 9

    def test_one_prg_with_even_variants(self, op, re):
        # Normalises odd end to even marker at site end
        re.side_effect = [to_bytes([5, 2, 6, 3, 5])]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        assert result == [5, 2, 6, 3, 6]

    def test_one_prg_with_direct_deletion(self, op, re):
        re.side_effect = [to_bytes([1, 2, 5, 6, 4, 5, 1])]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        assert result == [1, 2, 5, 6, 4, 6, 1]

    def test_even_marker_never_seen_fails(self, op, re):
        re.side_effect = [to_bytes([1, 2, 6, 3, 6, 4, 6])]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        with pytest.raises(PRGAggregationError):
            result = get_aggregated_prgs(agg, file_names)

    def test_site_marker_reused_fails(self, op, re):
        re.side_effect = [to_bytes([1, 2, 5, 3, 5, 4, 5, 3, 6, 2, 5])]
        file_names = {"prg_1": Path("dummy")}
        agg = PRGAggregator()
        with pytest.raises(PRGAggregationError):
            result = get_aggregated_prgs(agg, file_names)

    def test_two_prgs_correct_aggregation(self, op, re):
        prgs = [[1, 2, 5, 4, 4, 6, 3, 6, 1], [1, 3, 5, 1, 6, 1, 1, 6]]
        re.side_effect = map(to_bytes, prgs)
        file_names = {"prg_1": Path("dummy"), "prg_2": Path("dummo")}
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        expected = prgs[0] + [1, 3, 7, 1, 8, 1, 1, 8]
        assert result == expected
        assert agg.next_allocated == 9

    def test_multiple_prgs_one_nested_correct_aggregation(self, op, re):
        prgs = [
            [1, 2, 5, 4, 4, 6, 3, 6, 1],
            [1, 3, 5, 1, 7, 1, 1, 8, 1, 8, 2, 6, 6],
            [1, 2, 3, 4],
            [29, 2, 30, 4, 29],  # can use any var number
        ]
        re.side_effect = map(to_bytes, prgs)
        file_names = {
            "prg_1": Path("dummy"),
            "prg_2": Path("dummo"),
            "prg_3": Path("foo"),
            "prg_4": Path("foo"),
        }
        agg = PRGAggregator()
        result = get_aggregated_prgs(agg, file_names)
        expected = (
            prgs[0]
            + [1, 3, 7, 1, 9, 1, 1, 10, 1, 10, 2, 8, 8]
            + prgs[2]
            + [11, 2, 12, 4, 12]
        )
        assert result == expected
        assert agg.next_allocated == 13
