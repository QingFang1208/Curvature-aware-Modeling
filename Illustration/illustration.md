按钮说明：

<img src="image/menu.png" align="center">

操作说明：

1. <font color=#0099ff>对模型建立分割</font>

   在label模式下通过选点按钮在模型表面有顺序的选择一系列路径点，算法自动建立它们之间的最短路径

   <img src="image/selectPoints0.png" align="center">

   切回label模式可以进行旋转（按住鼠标左键），点击选点按钮添加新的路径点（最后一个点与最初的点重合）形成闭环

   <img src="image/selectPoints1.png" align="center">

   点击建立闭环分割的按钮

   <img src="image/segmentation.png" align="center">

2. <font color=#0099ff>对分割进行变形</font>

   选择变形类型，默认球形参数化

   <img src="image/deformKind.png" align="center">

   选择头部part

   <img src="image/selectPart.png" align="center">

   run 执行(点击refresh可以撤销这一次的变形)

   <img src="image/sphere.png" align="center">

assemble拼接(拼接之后目前的代码无法通过refresh回上一步状态)

<img src="image/assemble.png" align="center">