按钮说明：

![menu](.\image\menu.png)

操作说明：

1. <font color=#0099ff>对模型建立分割</font>

   在label模式下通过选点按钮在模型表面有顺序的选择一系列路径点，算法自动建立它们之间的最短路径

   ![menu](.\image\selectPoints0.png)

   切回label模式可以进行旋转（按住鼠标左键），点击选点按钮添加新的路径点（最后一个点与最初的点重合）形成闭环

   ![menu](.\image\selectPoints1.png)

   点击建立闭环分割的按钮

   ![menu](.\image\segmentation.png)

2. <font color=#0099ff>对分割进行变形</font>

   选择变形类型，默认球形参数化

   ![menu](.\image\deformKind.png)

   选择头部part

   ![menu](.\image\selectPart.png)

   run 执行(点击refresh可以撤销这一次的变形)

   ![menu](.\image\sphere.png)

assemble拼接(拼接之后目前的代码无法通过refresh回上一步状态)

![menu](.\image\assemble.png)