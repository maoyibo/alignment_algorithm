import sys
from typing import Dict, List, Set, Tuple, Union
from numpy.typing import NDArray

BASE_REPRESENTED: Dict[str, Set[str]] = {'A': {'A'}, 'a': {'a'},
                                         'C': {'C'}, 'c': {'c'},
                                         'G': {'G'}, 'g': {'g'},
                                         'T': {'T', 'U'}, 't': {'t', 'u'},
                                         'U': {'T', 'U'}, 'u': {'t', 'u'},
                                         'W': {'A', 'T', 'U'}, 'w': {'a', 't', 'u'},
                                         'S': {'C', 'G'}, 's': {'c', 'g'},
                                         'M': {'A', 'C'}, 'm': {'a', 'c'},
                                         'K': {'G', 'T', 'U'}, 'k': {'g', 't', 'u'},
                                         'R': {'A', 'G'}, 'r': {'a', 'g'},
                                         'Y': {'C', 'T', 'U'}, 'y': {'c', 't', 'u'},
                                         'B': {'C', 'G', 'T', 'U'}, 'b': {'c', 'g', 't', 'u'},
                                         'D': {'A', 'G', 'T', 'U'}, 'd': {'a', 'g', 't', 'u'},
                                         'H': {'A', 'C', 'T', 'U'}, 'h': {'a', 'c', 't', 'u'},
                                         'V': {'A', 'C', 'G'}, 'v': {'a', 'c', 'g'},
                                         'N': {'A', 'C', 'G', 'T', 'U'}, 'n': {'a', 'c', 'g', 't', 'u'},
                                         '-': set()}


def base_match(chr1: str, chr2: str) -> bool:
    chr1_repr = BASE_REPRESENTED[chr1.upper()]
    chr2_repr = BASE_REPRESENTED[chr2.upper()]
    return len(chr1_repr & chr2_repr) != 0


def print_matrix(matrix: Union[NDArray, List[List[int]]], seq1: str, seq2: str):
    print(
        ''.join([f"{str(x) if x!=sys.maxsize and x!=-sys.maxsize else '∞':>6s}" for x in '--'+seq2]))
    for (base, row) in zip('-'+seq1, matrix):
        print(''.join(
            f"{str(x) if x!=sys.maxsize and x!=-sys.maxsize else '∞':>6s}" for x in [base]+list(row)))


def align_report(path: List[Tuple[int, int, str, str, str]],
                       score: int):
    seq1_align = []
    seq2_align = []
    match_align = []
    for trace in path:
        seq1_align.append(trace[2])
        seq2_align.append(trace[3])
        match_align.append(trace[4])
    s1, s2 = path[0][:2]
    e1, e2 = path[-1][:2]
    print(f"{s1+1:<5d}{''.join(seq1_align)}{e1+1:>5d}\tScore: {score}")
    print(f"{'':5s}{''.join(match_align)}{'':>5s}")
    print(f"{s2+1:<5d}{''.join(seq2_align)}{e2+1:>5d}")


def align_report_gotoh(path: List[Tuple[str, int, int, str, str, str]],
                       score: int):
    seq1_align = []
    seq2_align = []
    match_align = []
    for trace in path:
        seq1_align.append(trace[3])
        seq2_align.append(trace[4])
        match_align.append(trace[5])
    s1, s2 = path[0][1:3]
    e1, e2 = path[-1][1:3]
    print(f"{s1+1:<5d}{''.join(seq1_align)}{e1+1:>5d}\tScore: {score}")
    print(f"{'':5s}{''.join(match_align)}{'':>5s}")
    print(f"{s2+1:<5d}{''.join(seq2_align)}{e2+1:>5d}")


def get_max_score_positions(matrix: Union[NDArray, List[List[int]]], max_score: int) -> List[Tuple[int, int]]:
    positions: List[Tuple[int, int]] = []
    for i, row in enumerate(matrix):
        for j, cell in enumerate(row):
            if cell == max_score:
                positions.append((i, j))
    return positions
