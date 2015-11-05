# ISAR image fusion

### Implement pedestrian detection and tracking by Hog feature and Kalman filter

### Detect by Hog feature
[Histograms of Oriented Gradients(HOG)](https://en.wikipedia.org/wiki/Histogram_of_oriented_gradients) feature extraction.
* Obtaining the human gradient information feature according to HoG operator.
* With overlapping block moving, the feature is a histograms of 4000 dimensions.

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
