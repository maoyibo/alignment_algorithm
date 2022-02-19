import sys
from typing import List, Tuple, Iterator
from alignment_algorithm.utils import base_match, get_max_score_positions, print_matrix, align_report_gotoh

import typer

app = typer.Typer()


def compute_scoring_matrix(seq1: str,
                           seq2: str,
                           match: int = 1,
                           mismatch: int = -1,
                           gap_open: int = -3,
                           gap_extend: int = -1) -> Tuple[List[List[int]], List[List[int]], List[List[int]], int]:
    scoring_matrix: List[List[int]] = []
    horizontal_matrix: List[List[int]] = []  # 水平Gap
    vertical_matrix: List[List[int]] = []  # 垂直Gap
    max_score = -sys.maxsize
    for i in range(len(seq1)+1):
        scoring_matrix.append([0]*(len(seq2)+1))
        horizontal_matrix.append([0]*(len(seq2)+1))
        vertical_matrix.append([0]*(len(seq2)+1))
        for j in range(len(seq2)+1):
            if i == 0:  # i 不变，j 增加 表示从水平的加减
                scoring_matrix[i][j] = 0  # 表示从序列开始的gap
                # 设置为最小值，表示不会选择该路径
                vertical_matrix[i][j] = -sys.maxsize if j != 0 else 0
                horizontal_matrix[i][j] = scoring_matrix[i][j]
            elif j == 0:  # j 不变， i 增加 表示垂直的加减
                scoring_matrix[i][j] = 0  # 表示从序列开始的gap
                # 设置为最小值，表示不会选择该路径
                horizontal_matrix[i][j] = -sys.maxsize if i != 0 else 0
                vertical_matrix[i][j] = scoring_matrix[i][j]
            else:
                if base_match(chr1=seq1[i-1],
                              chr2=seq2[j-1]) == True:
                    diagonal_value = scoring_matrix[i-1][j-1]+match
                else:
                    diagonal_value = scoring_matrix[i-1][j-1]+mismatch
                horizontal_matrix[i][j] = max(horizontal_matrix[i][j-1]+gap_extend,
                                              scoring_matrix[i][j-1]+gap_open+gap_extend)
                vertical_matrix[i][j] = max(vertical_matrix[i-1][j]+gap_extend,
                                            scoring_matrix[i-1][j]+gap_open+gap_extend)
                scoring_matrix[i][j] = max(diagonal_value,
                                           horizontal_matrix[i][j],
                                           vertical_matrix[i][j],
                                           0)
                if scoring_matrix[i][j] > max_score:
                    max_score = scoring_matrix[i][j]
    return scoring_matrix, horizontal_matrix, vertical_matrix, max_score


def get_horizontal_neighboured(scoring_matrix:  List[List[int]],
                               horizontal_matrix:  List[List[int]],
                               seq2: str,
                               i: int,
                               j: int,
                               gap_open: int = -3,
                               gap_extend: int = -1) -> List[Tuple[str, int, int, str, str, str]]:
    neighboured_list: List[Tuple[str, int, int, str, str, str]] = []
    if horizontal_matrix[i][j] == horizontal_matrix[i][j-1]+gap_extend:
        neighboured_list.append(('left', i, j-1, '-', seq2[j-1], ' '))
    if horizontal_matrix[i][j] == scoring_matrix[i][j-1]+gap_open+gap_extend:
        neighboured_list.append(('diagonal', i, j-1, '-', seq2[j-1], ' '))
    return neighboured_list


def get_vertical_neighboured(scoring_matrix:  List[List[int]],
                             vertical_matrix:  List[List[int]],
                             seq1: str,
                             i: int,
                             j: int,
                             gap_open: int = -3,
                             gap_extend: int = -1) -> List[Tuple[str, int, int, str, str, str]]:
    neighboured_list: List[Tuple[str, int, int, str, str, str]] = []
    if vertical_matrix[i][j] == vertical_matrix[i-1][j]+gap_extend:
        neighboured_list.append(('up', i-1, j, seq1[i-1], '-', ' '))
    if vertical_matrix[i][j] == scoring_matrix[i-1][j]+gap_open+gap_extend:
        neighboured_list.append(('diagonal', i-1, j, seq1[i-1], '-', ' '))
    return neighboured_list


def get_neighboured(scoring_matrix:  List[List[int]],
                    horizontal_matrix:  List[List[int]],
                    vertical_matrix:  List[List[int]],
                    seq1: str,
                    seq2: str,
                    matrix_type: str,  # diagonal, up, left
                    i: int,
                    j: int,
                    match: int = 1,
                    mismatch: int = -1,
                    gap_open: int = -3,
                    gap_extend: int = -1) -> List[Tuple[str, int, int, str, str, str]]:
    if matrix_type == 'up':  # i-1
        neighboured_list = get_vertical_neighboured(scoring_matrix=scoring_matrix,
                                                    vertical_matrix=vertical_matrix,
                                                    seq1=seq1,
                                                    i=i,
                                                    j=j,
                                                    gap_open=gap_open,
                                                    gap_extend=gap_extend)
    elif matrix_type == 'left':  # j-1:
        neighboured_list = get_horizontal_neighboured(scoring_matrix=scoring_matrix,
                                                      horizontal_matrix=horizontal_matrix,
                                                      seq2=seq2,
                                                      i=i,
                                                      j=j,
                                                      gap_open=gap_open,
                                                      gap_extend=gap_extend)
    else:
        assert matrix_type == 'diagonal', f"Unexpected matrix type: {matrix_type}"
        neighboured_list = []
        if i > 0 and j > 0:
            if base_match(chr1=seq1[i-1], chr2=seq2[j-1]) == True and scoring_matrix[i][j] == scoring_matrix[i-1][j-1] + match:
                neighboured_list.append(
                    ('diagonal', i-1, j-1, seq1[i-1], seq2[j-1], '|'))
            elif base_match(chr1=seq1[i-1], chr2=seq2[j-1]) == False and scoring_matrix[i][j] == scoring_matrix[i-1][j-1] + mismatch:
                neighboured_list.append(
                    ('diagonal', i-1, j-1, seq1[i-1], seq2[j-1], '*'))
        # i-1 up, seq2 gap
        if i > 0 and scoring_matrix[i][j] == vertical_matrix[i][j]:
            neighboured_list.extend(get_vertical_neighboured(scoring_matrix=scoring_matrix,
                                                             vertical_matrix=vertical_matrix,
                                                             seq1=seq1,
                                                             i=i,
                                                             j=j,
                                                             gap_open=gap_open,
                                                             gap_extend=gap_extend))
        # j-1 left, seq1 gap
        if j > 0 and scoring_matrix[i][j] == horizontal_matrix[i][j]:
            neighboured_list.extend(get_horizontal_neighboured(scoring_matrix=scoring_matrix,
                                                               horizontal_matrix=horizontal_matrix,
                                                               seq2=seq2,
                                                               i=i,
                                                               j=j,
                                                               gap_open=gap_open,
                                                               gap_extend=gap_extend))
    return neighboured_list


def local_traceback(scoring_matrix:  List[List[int]],
                    horizontal_matrix:  List[List[int]],
                    vertical_matrix:  List[List[int]],
                    seq1: str,
                    seq2: str,
                    matrix_type: str,
                    i: int,
                    j: int,
                    match: int = 1,
                    mismatch: int = -1,
                    gap_open: int = -3,
                    gap_extend: int = -1) -> Iterator[List[Tuple[str, int, int, str, str, str]]]:
    if (i == 0 and j == 0) or scoring_matrix[i][j] == 0:
        yield []
    else:
        for neighboured in get_neighboured(scoring_matrix=scoring_matrix,
                                           horizontal_matrix=horizontal_matrix,
                                           vertical_matrix=vertical_matrix,
                                           seq1=seq1,
                                           seq2=seq2,
                                           matrix_type=matrix_type,
                                           i=i,
                                           j=j,
                                           match=match,
                                           mismatch=mismatch,
                                           gap_open=gap_open,
                                           gap_extend=gap_extend):
            for path in local_traceback(scoring_matrix=scoring_matrix,
                                        horizontal_matrix=horizontal_matrix,
                                        vertical_matrix=vertical_matrix,
                                        seq1=seq1,
                                        seq2=seq2,
                                        matrix_type=neighboured[0],
                                        i=neighboured[1],
                                        j=neighboured[2],
                                        match=match,
                                        mismatch=mismatch,
                                        gap_open=gap_open,
                                        gap_extend=gap_extend):
                yield path+[neighboured]


@app.command()
def gotoh_local(seq1: str,
                seq2: str,
                match: int = 1,
                mismatch: int = -1,
                gap_open: int = -3,
                gap_extend: int = -1,
                show_matrix: bool = False,
                show_report: bool = True,
                show_path: bool = False):
    """Trace function
    """
    scoring_matrix, horizontal_matrix, vertical_matrix, max_score = compute_scoring_matrix(seq1=seq1,
                                                                                           seq2=seq2,
                                                                                           match=match,
                                                                                           mismatch=mismatch,
                                                                                           gap_open=gap_open,
                                                                                           gap_extend=gap_extend)
    if show_matrix == True:
        print('Vertical matrix:')
        print_matrix(matrix=vertical_matrix,
                     seq1=seq1,
                     seq2=seq2)
        print('Scoring matrix:')
        print_matrix(matrix=scoring_matrix,
                     seq1=seq1,
                     seq2=seq2)
        print('Horizontal matrix:')
        print_matrix(matrix=horizontal_matrix,
                     seq1=seq1,
                     seq2=seq2)
    for (i, j) in get_max_score_positions(matrix=scoring_matrix,
                                          max_score=max_score):
        for path in local_traceback(scoring_matrix=scoring_matrix,
                                    horizontal_matrix=horizontal_matrix,
                                    vertical_matrix=vertical_matrix,
                                    seq1=seq1,
                                    seq2=seq2,
                                    matrix_type='diagonal',
                                    i=i,
                                    j=j,
                                    match=match,
                                    mismatch=mismatch,
                                    gap_open=gap_open,
                                    gap_extend=gap_extend):
            if show_path:
                print(path)
            if show_report:
                align_report_gotoh(path,
                                   score=max_score)


if __name__ == "__main__":
    app()
