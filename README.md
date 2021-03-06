# ISAR image fusion

### Implement ISAR image fusion from different observing angles.
Carry out range alignment on backscattered data, and then implement phase correction. After extracting scattering centers, estimate aspect angle and rotation rate by using slope function and triangle method. Further implement image fusion which can promote image resolution in actual operation, thus obtaining more information about targets.

### Scattering Center Model: the basis of radar imaging and Target Model are as below
<img src="/resofmid/target.png" alt="4" width="300px"/>
<img src="/resultimage/rangealign_phsecorrection/originalmodel.png" alt="4" width="300px"/>


### Small Bandwith and Small Angle ISAR Imaging Simulation(Rotation Center Model)
Simple Motion Compensation method: CrossCorrelation, before and afer are as below。

<img src="/resofmid/1.png" alt="4" width="300px"/>
<img src="/resofmid/3.png" alt="4" width="300px"/>

### Big Bandwidth and Big Angle
Deal with electromagnetic data: PFA method , before and after are as below

<img src="/resofmid/4.png" alt="4" width="300px"/>
<img src="/resofmid/44.png" alt="4" width="300px"/>

### Range Alignment
Carry out range motion estimation and range alignment with CrossCorrelation method :
Range Profile of target Before and After are as below

<img src="/resofmid/rangepro1.png" alt="4" width="300px"/>
<img src="/resofmid/rangepro2.png" alt="4" width="300px"/>

### Rotational velocity estimationss
* Radon Detection of Lines method, Azimuth Rescaling before is as below, from which we can obtain the 
correct size of target

<img src="/resofmid/7.png" alt="7" width="300px"/>


* Simulated azimuth rotation is 0.04rad/s, and estimated velocity is 0.0360rad/s, and the rescaling
result is as below

<img src="/resofmid/8.png" alt="8" width="300px"/>

### Fusion example
* After rotaional velocity estimation and azimuth angle estimation , axis transformation , the result
of Fusion is as below, as can be seen , problem exists which needs to be solved

<img src="/resofmid/fusion_example.png" alt="fusion_example" width="700px"/>
