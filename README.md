# Curvature-aware-Modeling
This is a C++ implementation of the curvature modeling tool in the following paper:

Qing Fang, Zheng-Yu Zhao, Zhong-Yuan Liu, [Ligang Liu](http://staff.ustc.edu.cn/~lgliu/), [Xiao-Ming Fu](http://staff.ustc.edu.cn/~fuxm/). [Metric First Reconstruction for Interactive Curvature-aware Modeling](https://rec.ustc.edu.cn/share/ca848e70-82ce-11ea-b86e-171a380b3613). *Computer-Aided Design (SPM)*, 2020.

The code is written by Qing Fang and Zheng-Yu Zhao using Microsoft Visual Studio 2017 based on the [Surface Mesh Framework](http://staff.ustc.edu.cn/~fuxm/code/index.html#sec_surface_framework).  It contains:

1.  A rough segmentation tool to partition mesh into different parts.
2.  Three modeling operations for selected parts.

## External Libraries

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (3.2)

- [SuitSparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows)

## Data Format

*.seg :  face label array. The number records which part the face belongs.
