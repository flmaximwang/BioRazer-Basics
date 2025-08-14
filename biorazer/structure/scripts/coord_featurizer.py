import numpy as np
from sklearn.decomposition import PCA

def pca_coords(coords):
    '''
    使用 PCA 计算一组三维坐标的方向向量, 可用于计算
    - 螺旋的方向
    
    向量的范数为 1, 方向从第一个坐标指向最后一个坐标.
    '''
    
    coords = np.array(coords)
    if not coords.shape[1] == 3:
        raise ValueError("Coordinates must be a 2D array with shape (n_samples, 3).")
    if coords.shape[0] < 2:
        raise ValueError("Not enough coordinates to calculate direction vector.")
    
    # 使用 PCA 计算方向向量
    pca = PCA()
    pca.fit(coords)
    direction_vector = pca.components_[0]
    direction_vector = direction_vector / np.linalg.norm(direction_vector)
    direction_vector = np.sign(np.dot(direction_vector, coords[-1] - coords[0])) * direction_vector
    
    return direction_vector