from pathlib import Path
from typing import List


class Deletion:
    """
    Interval-based representation of a deletion record.
    start and stop are both 1-based inclusive.
    """

    def __init__(self, start: int, stop: int, del_len: int, samples: set):
        self.start = start
        self.stop = stop
        self.del_len = del_len
        self.samples = samples

    def modify_by(self, left_delta: int, right_delta: int):
        self.start += left_delta
        self.stop += right_delta

    def _contains(self, position: int):
        return self.start <= position <= self.stop

    def spans(self, other: "Deletion"):
        return self._contains(other.start) and self._contains(other.stop)

    def overlaps(self, other: "Deletion"):
        return self._contains(other.start) or self._contains(other.stop)

    def __len__(self) -> int:
        return self.stop - self.start + 1

    def __lt__(self, other: "Deletion") -> bool:
        return self.start < other.start

    def __eq__(self, other: "Deletion") -> bool:
        return self.start == other.start and self.stop == other.stop

    def __repr__(self):
        return f"[{self.start}, {self.stop}]: {self.samples}"


Deletions = List[Deletion]


def load_input_dels(input_dels_bed) -> Deletions:
    """
    Loads Deletions from a bed file
    """
    input_dels: Deletions = list()
    with Path(input_dels_bed).open() as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            rows = line.split("\t")
            samples = set(rows[3].strip().split(","))
            start, stop = int(rows[1]) + 1, int(rows[2])  # +1 : bed start is 0-based
            del_len = stop - start  # That's how I picked them: len(alt) == 1
            input_dels.append(Deletion(start, stop, del_len, samples))

    return input_dels
