from io import StringIO

import pytest

from make_prgs.extend_intervals import (
    DisjointInterval,
    load_existing_features,
    extend_features,
)


class TestLoadFeatures:
    def test_IntervalEquality(self):
        first = DisjointInterval(1, 2, "ref1")
        second = DisjointInterval(1, 2, "ref1")
        assert first == second

    def test_load_two_overlapping_intervals_fails(self):
        feature_stream = "ref1\t1\t10\tgene1\n" "ref1\t3\t12\tgene2\n"
        with pytest.raises(ValueError):
            result = load_existing_features(StringIO(feature_stream))

    def test_load_two_intervals_two_features(self):
        feature_stream = "ref1\t1\t10\tgene1\n" "ref2\t1\t10\tgene1\n"
        result = load_existing_features(StringIO(feature_stream))
        expected = {
            "ref1": [DisjointInterval(1, 10, "gene1")],
            "ref2": [DisjointInterval(1, 10, "gene1")],
        }
        assert result == expected

    def test_load_intervals_in_same_features(self):
        feature_stream = (
            "ref1\t1\t10\tgene1\n" "ref1\t20\t40\tgene2\n" "ref2\t1\t10\tgene1\n"
        )
        result = load_existing_features(StringIO(feature_stream))
        expected = {
            "ref1": [
                DisjointInterval(1, 10, "gene1"),
                DisjointInterval(20, 40, "gene2"),
            ],
            "ref2": [DisjointInterval(1, 10, "gene1")],
        }
        assert result == expected


@pytest.fixture
def feature_pool():
    pool = {
        "ref1": [
            DisjointInterval(1, 10, "gene1"),
            DisjointInterval(20, 40, "gene2"),
            DisjointInterval(200, 300, "gene3"),
        ],
        "ref2": [DisjointInterval(1, 100, "gene4")],
    }
    return pool


class TestExtendFeatures:
    def test_extend_single_feature(self, feature_pool):
        feature_pool.pop("ref1")
        extend_features(feature_pool, 100)
        expected = {"ref2": [DisjointInterval(0, 200, "gene4")]}
        assert feature_pool == expected

    def test_extend_two_features_no_overlap(self, feature_pool):
        feature_pool.pop("ref2")
        extend_features(feature_pool, 3)
        expected = {
            "ref1": [
                DisjointInterval(0, 13, "gene1"),
                DisjointInterval(17, 43, "gene2"),
                DisjointInterval(197, 303, "gene3"),
            ]
        }
        assert feature_pool == expected

    def test_extend_features_with_overlap(self, feature_pool):
        feature_pool.pop("ref2")
        extend_features(feature_pool, 100)
        expected = {
            "ref1": [
                DisjointInterval(0, 20, "gene1"),
                DisjointInterval(20, 140, "gene2"),
                DisjointInterval(140, 400, "gene3"),
            ]
        }
        assert feature_pool == expected
