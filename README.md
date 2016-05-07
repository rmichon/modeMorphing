# Mode Morphing and 3d Printing

## Bibliography Review

* **Henri Penttinen and Matti Karjalainen and Aki Härmä, Morphing Instrument Body Models**: Kind of the inspiration for our project. The mode interpolation technique that we use is a little bit similar to the one they describe (but not quite). I need to completely re-read this paper to compare it to what we're doing so that we can properly reference it in our ASA paper...

* **Li, Dingzeyu and Levin, David I.W. and Matusik, Wojciech and Zheng, Changxi: Acoustic Voxels: Computational Optimization of Modular Acoustic Filters**: I haven't read it yet but seems very linked to what we're doing. Also, their reference section probably has a bunch of stuff worth investigating (TODO).

 * **Gaurav Bharaj and David Levin and James Tompkin and Yun Fei and Hanspeter Pfister and Wojciech Matusik and Changxi Zheng: Computational Design of Metallophone Contact Sounds**: I haven't read it but it looks very promising and related to our work. This thing: http://www.reed.edu/reed_magazine/sallyportal/posts/2016/physics-major-invents-new-kind-of-bell.html that Gina recently posted on local-users is kind of similar too. 

* **A potentially interesting reference: http://aleph1.audio/** about a new kind of 3d printed speaker and acoustic measurements for it.

* **Other potentially interesting references: http://www.roboticstomorrow.com/article/2015/11/3d-printing-and-acoustics-rapid-prototyping-of-sound-diffusers/7096, http://www.materialise.com/cases/the-peugeot-fractal-concept-car-3d-printing-acoustic-interiors, http://scitation.aip.org/content/aip/journal/apl/108/6/10.1063/1.4941338, http://www.hovalin.com/, http://3dprintingindustry.com/2014/12/21/xl4d-3d-printed-smartphone-case/**


## General Notes/Observations

### Method (Juste to Summarize)

* Create a 2d matrix of resonators. The boundaries of this matrix are: 
	* square resonators <-> round resonators 
	* small resonators <-> large resonators
* Design the CAD models of the 4 corners of this matrix (large square, small square, large round and small round: "set 1")
* Design the CAD models of all the intermediate states (medium square, medium round, large square/round, medium square/round and small square/round: "set 2")
* Print all of them
* Measure their impulse response
* Extract their frequency response from their impulse response
* Detect all the potential modes of all the objects of set 1 between 0 and 5KHz. For that, we considered any peak in the spectrum as a potential mode (no threshold). The reason for that is that a very significant mode on one corner of the matrix can potentially become totally obsolete on the other side but we still need to be able to track it.
* NOTE: for the following steps, we did things manually because we had to rush but in practice (now that we know the results), we know that this can be done automatically. Thus, when I'll use the word "automatically" for the following steps, it might mean that this is something that we did manually in practice... For the ASA paper, we should say that things were done automatically because it makes more sense scintifically speaking :)... 
* 

### Others

* Don't forget to talk about the large data set and machine learning idea.

## Content

TODO
