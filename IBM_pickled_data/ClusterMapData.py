from typing import List

import numpy as np


class ClusterMapData(object):

    def __init__(self,
                 data: np.array,
                 row_linkage: np.ndarray,
                 col_linkage: np.ndarray,
                 row_labels: List[str],
                 col_labels: List[str],
                 col_colors: List[str],
                 col_color_indices: List[int]) -> None:
        self.data = data
        self.row_linkage = row_linkage
        self.col_linkage = col_linkage
        self.row_labels = row_labels
        self.col_labels = col_labels
        self.col_colors = col_colors
        self.col_color_indices = col_color_indices
