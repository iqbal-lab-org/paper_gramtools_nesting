from pathlib import Path
from typing import List


class Interval:
    """
    start and stop are both 1-based inclusive.
    """

    def __init__(self, start: int, stop: int, samples: set = None, size: int = None):
        self.start = start
        self.stop = stop
        self.samples = samples

        self.size = size
        self.size = self.__len__()

    def modify_by(self, left_delta: int, right_delta: int):
        self.start += left_delta
        self.stop += right_delta

    def _contains(self, position: int):
        return self.start <= position <= self.stop

    def spans(self, other: "Interval"):
        return self._contains(other.start) and self._contains(other.stop)

    def overlaps(self, other: "Interval"):
        return self._contains(other.start) or self._contains(other.stop)

    def __len__(self) -> int:
        if self.size is None:
            self.size = self.stop - self.start + 1
        return self.size

    def __lt__(self, other: "Interval") -> bool:
        return self.start < other.start

    def __eq__(self, other: "Interval") -> bool:
        return self.start == other.start and self.stop == other.stop

    def __repr__(self):
        return f"[{self.start}, {self.stop}]: {self.samples}"


Intervals = List[Interval]


def load_input_dels(input_dels_bed) -> Intervals:
    """
    Loads deletions as Intervals from a bed file
    """
    input_dels: Intervals = list()
    with Path(input_dels_bed).open() as f:
        for i, line in enumerate(f):
            rows = line.split("\t")
            samples = set(rows[3].strip().split(","))
            start, stop = int(rows[1]) + 1, int(rows[2])  # +1 : bed start is 0-based
            del_len = stop - start  # That's how I picked them: len(alt) == 1
            input_dels.append(Interval(start, stop, samples, del_len))

    return input_dels
