QSAP Instances used in: 

- Silva, Allyson, Leandro C. Coelho, and Maryam Darvish. 
"Quadratic assignment problem variants: A survey and an effective parallel memetic iterated tabu search." 
European Journal of Operational Research 292.3 (2021): 1066-1084.

- Cordeau, Jean-François, et al. 
"A memetic heuristic for the generalized quadratic assignment problem." 
INFORMS Journal on Computing 18.4 (2006): 433-443.

Each instance folder contains four files Theta,T,A and D.

* The theta file corresponds to the metadata of the instance, the Number of sites n, the number of facilities m and finally the traffic unit cost w
* T corresponds to the traffic intensity between facilities, and it is a (n,n) matrix showed in long format. 
  The file contains three columns the first two are the indexes of the facilities, and the third column corresponds
  to the traffic intesity value between the facilities.
* D correponds to the distance matrix between sites, and it is a (n,n) matrix showed in long format.
  The file contains three columns the first two are the indexes of the sites, and the third column corresponds to the distance between sites.
* A corresoonds to the cost of assigning a facility to a site and it is a (n,m) matrix showed in long format.
  The file contains three columns the first two are the indexes of the facilities and sites, and the third column corresponds to the 
  cost of assignment of facilities to sites.

Please note that all indexes start counting from 1 and not 0
