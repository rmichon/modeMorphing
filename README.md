# Mode Morphing and 3d Printing

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

## Content

TODO
