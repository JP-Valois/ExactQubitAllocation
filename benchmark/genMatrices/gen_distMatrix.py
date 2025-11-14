from utils import * # type: ignore
import numpy as np


def BFS_distances (Graph:list[list], start:int) -> np.array:

    m = len(Graph)

    dist = np.full(m, np.inf)

    dist[start] = -1
    queue = [start]

    while (len(queue) > 0):

        node = queue.pop(0)

        for neighbor in Graph[node]:

            if dist[neighbor] == np.inf:
                dist[neighbor] = dist[node] + 1
                queue.append(neighbor)
            
    dist[start] = 0

    return dist


def DistanceMatrix (Graph:list[list]) -> np.array:

    m = len(Graph)

    D = np.zeros((m,m), dtype=int)

    for i in range(m):
        D[i] = BFS_distances(Graph, i)

    return D


dir_name = "dist"


## 16-qubit IBM Melbourne NISQ
Graph_Melbourne16 = [[1, 15],
                     [0, 14, 2],
                     [1, 13, 3],
                     [2, 12, 4],
                     [3, 11, 5],
                     [4, 10, 6],
                     [5, 9, 7],
                     [6, 8],
                     [9, 7],
                     [8, 6, 10],
                     [9, 5, 11],
                     [10, 4, 12],
                     [11, 3, 13],
                     [12, 2, 14],
                     [13, 1, 15],
                     [0, 14]]


output_file_path_1 = "../" + dir_name + "/16_melbourne.csv"
SaveMatrix_csv(DistanceMatrix(Graph_Melbourne16), output_file_path_1)
print("Saved distance matrix as", output_file_path_1)


## 20-qubit IBM Boeblingen NISQ
Graph_Boeblingen20 = [[1],           # 0
                      [0, 2, 6],     # 1
                      [1, 3],        # 2
                      [2, 4, 8],     # 3
                      [3],           # 4
                      [6, 10],       # 5
                      [1, 5, 7],     # 6
                      [6, 8, 12],        # 7
                      [3, 7, 9],     # 8
                      [8, 14],       # 9
                      [5, 11],       # 10
                      [10, 12, 16],  # 11
                      [7, 11, 13],   # 12
                      [12, 14, 18],  # 13
                      [9, 13],       # 14
                      [16],          # 15
                      [11, 15, 17],  # 16
                      [16, 18],      # 17
                      [13, 17, 19],  # 18
                      [18]]          # 19


output_file_path_2 = "../" + dir_name + "/20_boeblingen.csv"
SaveMatrix_csv(DistanceMatrix(Graph_Boeblingen20), output_file_path_2)
print("Saved distance matrix as", output_file_path_2)


## 27-qubit IBM Cairo NISQ
Graph_Cairo27 = [[1],           # 0
                 [0, 2, 4],     # 1
                 [1, 3],        # 2
                 [2, 5],        # 3
                 [1, 7],        # 4
                 [3, 8],        # 5
                 [7],           # 6
                 [4, 6, 10],    # 7
                 [5, 9, 11],    # 8
                 [8],           # 9
                 [7, 12],       # 10
                 [8, 14],       # 11
                 [10, 13, 15],  # 12
                 [12, 14],      # 13
                 [11, 13, 16],  # 14
                 [12, 18],      # 15
                 [14, 19],      # 16
                 [18],          # 17
                 [15, 17, 21],  # 18
                 [16, 20, 22],  # 19
                 [19],          # 20
                 [18, 23],      # 21
                 [19, 25],      # 22
                 [21, 24],      # 23
                 [23, 25],      # 24
                 [22, 24, 26],  # 25
                 [25]]          # 26


output_file_path_3 = "../" + dir_name + "/27_cairo.csv"
SaveMatrix_csv(DistanceMatrix(Graph_Cairo27), output_file_path_3)
print("Saved distance matrix as", output_file_path_3)
