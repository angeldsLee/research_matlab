# ISAR image fusion

### Implement ISAR image fusion from different observing angles.
Carry out range alignment on backscattered data, and then implement phase correction. After extracting scattering centers, estimate aspect angle and rotation rate by using slope function and triangle method. Further implement image fusion which can promote image resolution in actual operation, thus obtaining more information about targets.

### Scatter model: the basis of radar imaging
<img src="/resultimage/rangealign_phsecorrection/originalmodel.png" alt="Hog1" width="300px"/>

<img src="/features/capture1.png" alt="Hog1" width="300px"/>
<img src="/features/capture2.png" alt="Hog2" width="300px"/>
<img src="/features/capture3.png" alt="Hog3" width="300px"/>
<img src="/features/capture4.png" alt="Hog4" width="300px"/>

* The extracted features are use by [SVM]() to detection pedestrian candidates.

### Tracking by Kalman filter
* Using [Kalman filter](https://en.wikipedia.org/wiki/Kalman_filter) to predict the movement of pedestrian.
* The movement vectors look like as following figures.

<img src="/figures/capture8.png" alt="move1" width="300px"/>
<img src="/figures/5555.png" alt="move3" width="300px"/>
<img src="/figures/capture10.png" alt="move3" width="300px"/>
<img src="/figures/capture9.png" alt="move2" width="300px"/>
