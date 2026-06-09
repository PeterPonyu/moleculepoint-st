from __future__ import annotations
from math import hypot
from .contracts import MoleculePoint, PointPatternSummary


def nearest_neighbor_distance(points: list[MoleculePoint], gene: str) -> float:
    selected = [p for p in points if p.gene == gene]
    if len(selected) < 2:
        return float("inf")
    vals = []
    for i, p in enumerate(selected):
        vals.append(min(hypot(p.x-q.x, p.y-q.y) for j, q in enumerate(selected) if i != j))
    return sum(vals) / len(vals)


def summarize_gene(points: list[MoleculePoint], gene: str, low_count_floor: int = 5) -> PointPatternSummary:
    count = sum(1 for p in points if p.gene == gene)
    return PointPatternSummary(gene, count, nearest_neighbor_distance(points, gene), count < low_count_floor)


def colocalization_fraction(points: list[MoleculePoint], gene_a: str, gene_b: str, radius: float = 0.35) -> float:
    a = [p for p in points if p.gene == gene_a]
    b = [p for p in points if p.gene == gene_b]
    if not a or not b:
        return 0.0
    hits = sum(any(hypot(p.x-q.x, p.y-q.y) <= radius for q in b) for p in a)
    return hits / len(a)
