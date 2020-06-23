from unittest.mock import MagicMock

import pytest
from igraph import Graph

from msps_dimorphism.get_site_diversity_graphs import (
    Region,
    is_in_region,
    get_sample_indices,
    wire,
    make_site_graph,
)


@pytest.fixture(scope="class")
def site_json_data():
    return {"SEG": "seg1", "POS": "10"}


class TestRegion:
    def test_is_region_on_empty_region_returns_true(self, site_json_data):
        assert is_in_region(site_json_data, Region())

    def test_not_in_region(self, site_json_data):
        for region in [
            Region("seg2", 10, 20),
            Region("seg1", 11, 20),
            Region("seg1", 9, 9),
        ]:
            assert not is_in_region(site_json_data, region)

    def test_in_region(self, site_json_data):
        for region in [Region("seg1", 10, 10), Region("seg1", 8, 20)]:
            assert is_in_region(site_json_data, region)


class TestSampleIndices:
    jvcf = {"Samples": [{"Name": "s1"}, {"Name": "s2"}, {"Name": "s3"}]}

    def test_get_all_samples(self):
        result = get_sample_indices(self.jvcf)
        assert result == [0, 1, 2]

    def test_get_subset_of_samples(self):
        result = get_sample_indices(self.jvcf, ["s1", "s3"])
        assert result == [0, 2]


@pytest.fixture(scope="function")
def test_wire_data(request):
    jvcf = {
        "Child_Map": {"0": {"0": [1, 2], "1": [3]}, "1": {"2": [4]}},
        "Sites": MagicMock(),
    }

    graph = Graph(directed=True)
    num_sites = 6
    graph.add_vertices(num_sites)
    graph.vs["POS"] = [0] * num_sites
    graph.vs["populated"] = [False] * num_sites
    graph.vs["nesting_lvl"] = [1] * num_sites

    request.cls.jvcf = jvcf
    request.cls.graph = graph


class TestWire:
    def test_wire_no_nesting(self, test_wire_data):
        wire(2, 5, 1, self.jvcf, self.graph)
        assert sum(self.graph.vs["populated"]) == 1
        assert len(self.graph.es) == 1

    def test_wire_with_nesting(self, test_wire_data):
        wire(0, 5, 1, self.jvcf, self.graph)
        assert sum(self.graph.vs["populated"]) == 5
        assert self.graph.vs["nesting_lvl"] == [1, 2, 2, 2, 3, 1]
        assert self.graph.degree(mode="in") == [0, 1, 1, 1, 1, 2]
        assert self.graph.degree(mode="out") == [2, 1, 1, 1, 1, 0]


@pytest.fixture(scope="function")
def site_graph_data(request):
    num_sites = 8
    jvcf = {
        "Child_Map": {"0": {"0": [1], "1": [2]}, "4": {"1": [5]}, "5": {"1": [6]}},
        "Sites": [],
        "Lvl1_Sites": [0, 3, 4, 7],
    }
    for i in range(num_sites):
        jvcf["Sites"].append({"SEG": "seg1", "POS": i + 1})

    graph = Graph(directed=True)
    graph.add_vertices(num_sites)
    graph.vs["POS"] = [0] * num_sites
    graph.vs["populated"] = [False] * num_sites
    graph.vs["nesting_lvl"] = [1] * num_sites

    request.cls.jvcf = jvcf
    request.cls.graph = graph


class TestMakeSiteGraph:
    def test_one_site_only(self, site_graph_data):
        self.jvcf["Sites"] = [self.jvcf["Sites"][0]]
        self.jvcf["Child_Map"] = {}
        self.jvcf["Lvl1_Sites"] = [0]

        result = make_site_graph(self.jvcf, Region())

        assert len(result.vs) == 1

    def test_no_sites_in_region(self, site_graph_data):
        region = Region("invalid_seg", 1, 10)
        with pytest.raises(IndexError):
            result = make_site_graph(self.jvcf, region)

    def test_apply_region_filter(self, site_graph_data):
        region = Region("seg1", 4, 4)
        result = make_site_graph(self.jvcf, region)
        assert sum(result.vs["populated"]) == 1

    def test_apply_region_filter2(self, site_graph_data):
        # Includes a site, at index 4 and POS 5,
        # which has nested sites below it; these get processed too.
        region = Region("seg1", 4, 5)
        result = make_site_graph(self.jvcf, region)
        assert sum(result.vs["populated"]) == 4
        assert result.vs["nesting_lvl"] == [0, 0, 0, 1, 1, 2, 3, 0]
        assert result.degree(mode="in") == [0, 0, 0, 0, 1, 1, 1, 1]
        assert result.degree(mode="out") == [0, 0, 0, 1, 1, 1, 1, 0]

    def test_make_whole_graph(self, site_graph_data):
        result = make_site_graph(self.jvcf, Region())
        assert result.vs["nesting_lvl"] == [1, 2, 2, 1, 1, 2, 3, 1]
        assert result.degree(mode="in") == [0, 1, 1, 2, 1, 1, 1, 1]
        assert result.degree(mode="out") == [2, 1, 1, 1, 1, 1, 1, 0]
