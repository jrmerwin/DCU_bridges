# Extracted unchanged from the original notebook import and structure cells.
from collections import Counter
from itertools import combinations
from dataclasses import dataclass, asdict
import random
import sys

print("Python:", sys.version.split()[0])
print("DCU Rosetta notebook: native counting model, revised manuscript 0.5")

class DCUStructure:
    """Unordered parent pairs, exact ancestry, and the paper's recording cost."""

    def __init__(self):
        # IDs are storage labels only: 0 = a, 1 = b.
        self.parents = [None, None]
        self.chains = [1, 1]
        self.ancestors = [frozenset({0}), frozenset({1})]
        self.recorded = set()
        self.pair_to_id = {}

    def __len__(self):
        return len(self.parents)

    def path_weight(self, node):
        return 2 * self.chains[node] - 2

    def _checked_pairs(self, pairs):
        result = []

        for u, v in pairs:
            if type(u) is not int or type(v) is not int:
                raise TypeError("Parents must be integer IDs in this structure.")

            if not (0 <= u < len(self) and 0 <= v < len(self)):
                raise ValueError(
                    "Both parents must already exist before the burst."
                )

            if u == v:
                raise ValueError(
                    "A parent pair must contain two distinct objects."
                )

            pair = (min(u, v), max(u, v))

            if pair in self.pair_to_id:
                raise ValueError(
                    "This pair already exists; it cannot create a duplicate."
                )

            result.append(pair)

        if len(result) != len(set(result)):
            raise ValueError("A burst cannot contain the same pair twice.")

        return sorted(result)

    def _cost(self, pairs):
        hits = Counter()

        for u, v in pairs:
            # Shared ancestry is counted ONCE within each event.
            hits.update(self.ancestors[u] | self.ancestors[v])

        repeat_part = sum(
            2 * self.path_weight(z) * count
            for z, count in hits.items()
        )

        first_use_extra = sum(
            9 * self.path_weight(z)
            for z in hits
            if z not in self.recorded
        )

        return repeat_part + first_use_extra

    def recording_cost(self, pairs):
        """Quote the cost without changing the structure."""
        return self._cost(self._checked_pairs(pairs))

    def add_batch(self, pairs):
        """Commit a specified batch. All parents belong to the OLD structure."""
        pairs = self._checked_pairs(pairs)
        cost = self._cost(pairs)
        newborns = []

        for u, v in pairs:
            node = len(self)

            self.parents.append((u, v))
            self.chains.append(self.chains[u] + self.chains[v])

            self.ancestors.append(
                self.ancestors[u]
                | self.ancestors[v]
                | frozenset({node})
            )

            self.pair_to_id[(u, v)] = node
            self.recorded.update((u, v))
            newborns.append(node)

        return newborns, cost
