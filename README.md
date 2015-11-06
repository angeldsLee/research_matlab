# ISAR image fusion

### Implement ISAR image fusion from different observing angles.
Carry out range alignment on backscattered data, and then implement phase correction. After extracting scattering centers, estimate aspect angle and rotation rate by using slope function and triangle method. Further implement image fusion which can promote image resolution in actual operation, thus obtaining more information about targets.

### Scatter model: the basis of radar imaging
<img src="/resultimage/rangealign_phsecorrection/originalmodel.png" alt="scatter" width="300px"/>


### 小相干积累角度ISAR成像仿真(转台模型)
简单的运动补偿方法: 互相关法平动估计与补偿, 补偿前成像结果如下图所示。

<img src="/resofmid/1.png" alt="1" width="300px"/>

### 斜率估计和大相干积累角度PFA算法
仿真处理的实际数据:补偿前结果如下图所示。

<img src="/resofmid/4.png" alt="4" width="300px"/>

### Rotational velocity estimationss
* Radon Detection of Lines方法,方位向定标前如下图所示:

<img src="/resofmid/7.png" alt="7" width="300px"/>

* 实验假设转速为0.04rad/s,估计的转速为0.0360rad/s,从估计结果可见,基本反映了真是目标的尺寸信息。
方位向定标后如下图所示:

<img src="/resofmid/8.png" alt="8" width="300px"/>

### Fusion example
<img src="/resofmid/fusion_example.png" alt="fusion_example" width="700px"/>
