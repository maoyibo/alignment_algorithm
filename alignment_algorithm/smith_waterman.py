import sys
from typing import List, Tuple, Iterator
from alignment_algorithm.utils import base_match, get_max_score_positions, align_report, print_matrix

import typer

app = typer.Typer()


def compute_scoring_matrix(seq1: str,
                           seq2: str,
                           match: int = 1,
                           mismatch: int = -1,
                           gap_penalty: int = -2) -> Tuple[List[List[int]], int]:
    scoring_matrix: List[List[int]] = []
    max_score: int = -sys.maxsize
    for i in range(len(seq1)+1):
        scoring_matrix.append([0]*(len(seq2)+1))
        for j in range(len(seq2)+1):
            if i == 0:
                scoring_matrix[i][j] = 0
            elif j == 0:
                scoring_matrix[i][j] = 0
            else:
                if base_match(chr1=seq1[i-1],
                              chr2=seq2[j-1]) == True:
                    diagonal_value = scoring_matrix[i-1][j-1]+match
                else:
                    diagonal_value = scoring_matrix[i-1][j-1]+mismatch
                left_value = scoring_matrix[i][j-1]+gap_penalty
                up_value = scoring_matrix[i-1][j]+gap_penalty
                scoring_matrix[i][j] = max(diagonal_value,
                                           left_value,
                                           up_value,
                                           0)
                if scoring_matrix[i][j] > max_score:
                    max_score = scoring_matrix[i][j]
    return scoring_matrix, max_score


def get_neighboured(scoring_matrix:  List[List[int]],
                    seq1: str,
                    seq2: str,
                    i: int,
                    j: int,
                    match: int = 1,
                    mismatch: int = -1,
                    gap_penalty: int = -2) -> List[Tuple[int, int, str, str, str]]:
    neighboured_list: List[Tuple[int, int, str, str, str]] = []
    # diagonal
    if i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1] and scoring_matrix[i][j] == scoring_matrix[i-1][j-1] + match:
            neighboured_list.append((i-1, j-1, seq1[i-1], seq2[j-1], '|'))
        elif seq1[i-1] != seq2[j-1] and scoring_matrix[i][j] == scoring_matrix[i-1][j-1] + mismatch:
            neighboured_list.append((i-1, j-1, seq1[i-1], seq2[j-1], '*'))
    # i-1 up, seq2 gap
    if i > 0 and scoring_matrix[i][j] == scoring_matrix[i-1][j]+gap_penalty:
        neighboured_list.append((i-1, j, seq1[i-1], '-', ' '))
    # j-1 left, seq1 gap
    if j > 0 and scoring_matrix[i][j] == scoring_matrix[i][j-1]+gap_penalty:
        neighboured_list.append((i, j-1, '-', seq2[j-1], ' '))
    return neighboured_list


def local_traceback(scoring_matrix:  List[List[int]],
                    seq1: str,
                    seq2: str,
                    i: int,
                    j: int,
                    match: int = 1,
                    mismatch: int = -1,
                    gap_penalty: int = -2) -> Iterator[List[Tuple[int, int, str, str, str]]]:
    if (i == 0 and j == 0) or scoring_matrix[i][j] == 0:
        yield []
    else:
        for neighboured in get_neighboured(scoring_matrix=scoring_matrix,
                                           seq1=seq1,
                                           seq2=seq2,
                                           i=i,
                                           j=j,
                                           match=match,
                                           mismatch=mismatch,
                                           gap_penalty=gap_penalty):
            for path in local_traceback(scoring_matrix=scoring_matrix,
                                        seq1=seq1,
                                        seq2=seq2,
                                        i=neighboured[0],
                                        j=neighboured[1],
                                        match=match,
                                        mismatch=mismatch,
                                        gap_penalty=gap_penalty):
                yield path+[neighboured]


@app.command('tf')
def sw_trace_function(seq1: str,
                      seq2: str,
                      match: int = 1,
                      mismatch: int = -1,
                      gap_penalty: int = -2,
                      show_matrix: bool = False,
                      show_report: bool = True,
                      show_path: bool = False):
    """Trace function
    """
    scoring_matrix, max_score = compute_scoring_matrix(seq1=seq1,
                                                       seq2=seq2,
                                                       match=match,
                                                       mismatch=mismatch,
                                                       gap_penalty=gap_penalty)
    if show_matrix == True:
        print('Scoring matrix:')
        print_matrix(matrix=scoring_matrix,
                     seq1=seq1,
                     seq2=seq2)
    for (i, j) in get_max_score_positions(matrix=scoring_matrix,
                                          max_score=max_score):
        for path in local_traceback(scoring_matrix=scoring_matrix,
                                    seq1=seq1,
                                    seq2=seq2,
                                    i=i,
                                    j=j,
                                    match=match,
                                    mismatch=mismatch,
                                    gap_penalty=gap_penalty):
            if show_path:
                print(path)
            if show_report:
                align_report(path,
                             score=max_score)


def compute_scoring_trace_matrix(seq1: str,
                                 seq2: str,
                                 match: int = 1,
                                 mismatch: int = -1,
                                 gap_penalty: int = -2) -> Tuple[List[List[int]], List[List[int]], int]:
    # 1 for diagonal
    # 2 for up
    # 4 for left
    scoring_matrix: List[List[int]] = []
    trace_matrix: List[List[int]] = []
    max_score: int = 0
    for i in range(len(seq1)+1):
        scoring_matrix.append([0]*(len(seq2)+1))
        trace_matrix.append([0]*(len(seq2)+1))
        for j in range(len(seq2)+1):
            if i == 0:
                scoring_matrix[i][j] = 0
                trace_matrix[i][j] |= 4
            elif j == 0:
                scoring_matrix[i][j] = 0
                trace_matrix[i][j] |= 2
            else:
                if base_match(chr1=seq1[i-1],
                              chr2=seq2[j-1]) == True:
                    diagonal_value = scoring_matrix[i-1][j-1]+match
                else:
                    diagonal_value = scoring_matrix[i-1][j-1]+mismatch
                left_value = scoring_matrix[i][j-1]+gap_penalty
                up_value = scoring_matrix[i-1][j]+gap_penalty
                scoring_matrix[i][j] = max(diagonal_value,
                                           left_value,
                                           up_value,
                                           0)
                if left_value == scoring_matrix[i][j]:
                    trace_matrix[i][j] |= 4
                if up_value == scoring_matrix[i][j]:
                    trace_matrix[i][j] |= 2
                if diagonal_value == scoring_matrix[i][j]:
                    trace_matrix[i][j] |= 1
                if scoring_matrix[i][j] > max_score:
                    max_score = scoring_matrix[i][j]
    trace_matrix[0][0] = 0
    return scoring_matrix, trace_matrix, max_score


def global_traceback_trace_matrix(scoring_matrix:  List[List[int]],
                                  trace_matrix: List[List[int]],
                                  seq1: str,
                                  seq2: str,
                                  i: int,
                                  j: int) -> Iterator[List[Tuple[int, int, str, str, str]]]:
    if (i == 0 and j == 0) or scoring_matrix[i][j] == 0:
        yield []
    else:
        trace = trace_matrix[i][j]
        if trace & 1 == 1:  # diagonal
            for path in global_traceback_trace_matrix(scoring_matrix=scoring_matrix,
                                                      trace_matrix=trace_matrix,
                                                      seq1=seq1,
                                                      seq2=seq2,
                                                      i=i-1,
                                                      j=j-1):
                yield path+[(i-1, j-1, seq1[i-1], seq2[j-1], '|' if seq1[i-1] == seq2[j-1] else '*')]
        if trace & 2 == 2:  # up
            for path in global_traceback_trace_matrix(scoring_matrix=scoring_matrix,
                                                      trace_matrix=trace_matrix,
                                                      seq1=seq1,
                                                      seq2=seq2,
                                                      i=i-1,
                                                      j=j):
                yield path+[(i-1, j-1, seq1[i-1], '-', ' ')]
        if trace & 4 == 4:  # left
            for path in global_traceback_trace_matrix(scoring_matrix=scoring_matrix,
                                                      trace_matrix=trace_matrix,
                                                      seq1=seq1,
                                                      seq2=seq2,
                                                      i=i,
                                                      j=j-1):
                yield path+[(i-1, j-1, '-', seq2[j-1], ' ')]


@app.command('tm')
def sw_trace_matrix(seq1: str,
                    seq2: str,
                    match: int = 1,
                    mismatch: int = -1,
                    gap_penalty: int = -2,
                    show_matrix: bool = False,
                    show_report: bool = True,
                    show_path: bool = False):
    """Trace matrix
    """
    scoring_matrix, trace_matrix, max_score = compute_scoring_trace_matrix(seq1=seq1,
                                                                           seq2=seq2,
                                                                           match=match,
                                                                           mismatch=mismatch,
                                                                           gap_penalty=gap_penalty)
    if show_matrix == True:
        print('Scoring matrix:')
        print_matrix(matrix=scoring_matrix,
                     seq1=seq1,
                     seq2=seq2)
        print('Trace matrix:')
        print_matrix(matrix=trace_matrix,
                     seq1=seq1,
                     seq2=seq2)
    for (i, j) in get_max_score_positions(matrix=scoring_matrix,
                                          max_score=max_score):
        for path in global_traceback_trace_matrix(scoring_matrix=scoring_matrix,
                                                  trace_matrix=trace_matrix,
                                                  seq1=seq1,
                                                  seq2=seq2,
                                                  i=i,
                                                  j=j):
            if show_path:
                print(path)
            if show_report:
                align_report(path,
                             score=max_score)


if __name__ == "__main__":
    app()
