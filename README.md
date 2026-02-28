# law-minimum
Tests the Law of the Minimum in Agriculture

# TODO: 

* From Overleaf add quarto document in papers to create the versions of the theoretical graph
* Seperate script to pull county level NASS Data 
* Sepearte script to create county panel of soil measures and soil additives 
* Script for descriptive regressions of relationship of practices and soil measures 
* Script isolating GMM estimation fitting equation. 

Other Nick thoughts I have yet to address: 

* We need to figure out how to make the parameters conditional on 'soil health' (in our case OM)... your interaction term is a good approach. Let's go over these values more. 
* It might help our intuition if you make three versions of this graph
* One with three different alphas, one with three different betas, and one with three different gammas. 
* (a0 + a1*OM)(1 - (b0 + b1*OM)*exp(-(g0 + g1*OM))
* A caveat here is what does it mean for OM to be 0? This is not realistic. 
* Is this ok if our intercept is meaningless as long as the range of values over the domain of of observe OM makes sense?
* (a0 + a1*OM)(1 - (b0 + b1*OM)*P*exp(-(g0 + g1*OM)*P)