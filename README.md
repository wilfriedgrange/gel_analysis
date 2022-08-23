# gel_analysis

Automatically determines band intensities from an agarose gel's image (as a 8 bit file).

- For the demo, you may use the attached .txt file (below is the .jpg).

![IM1](https://user-images.githubusercontent.com/110093777/186207086-f66c5fc9-5ace-4b64-9f10-49bde83b6c26.jpg)

- You may have to change the following parameters for best results.

```
ws<-10 # smooth baseline (change !)
wm_lowMW<-10# wm for low Molecular Weight (change !)
wm_highMW<-15 #  wm for high Molecular Weight (change !)
half_size<-18 # half fitting range (change !)
m_cut<- 18 # used to separate lanes (change !)
```

- Here, the .R file should ouput 12 plots that correspond to the 12 lanes. 

For Lane 11

![11](https://user-images.githubusercontent.com/110093777/186207733-f9c7a6d8-27a1-400c-8e92-7f8a993f8ea1.jpeg)

For Lane 12

![12](https://user-images.githubusercontent.com/110093777/186207760-c2aa2b5e-6f64-48c3-babc-9dd1f4f1ef48.jpeg)
