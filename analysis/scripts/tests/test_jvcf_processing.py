import pytest

from jvcf_processing import Region, is_in_region, num_sites_under


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


@pytest.fixture(scope="class")
def child_map_data():
    return {"0": {"0": [1, 2], "1": [3, 4]}, "1": {"0": [5, 6, 7]}}


class TestSitesUnder:
    def test_nonparentsite_returns_zero(self, child_map_data):
        assert num_sites_under(child_map_data, "5") == 0

    def test_parentsite_returns_directchildren(self, child_map_data):
        assert num_sites_under(child_map_data, "1") == 3

    def test_parentsite_returns_allchildren_recursively(self, child_map_data):
        assert num_sites_under(child_map_data, "0") == 7
