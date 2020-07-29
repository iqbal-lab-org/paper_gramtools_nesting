import pytest

from msps_dimorphism.site_regions import Region, is_in_region


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
