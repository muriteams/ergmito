# ergmito 0.2-1-9999

* Added a `NEWS.md` file to track changes to the package.

* Fixed bug in `ergmito_boot`. The covariance matrix was estimated with the
  wrong sample, sometimes taking observations with missing.

* Better memory management: `ergmito_ptr` (a C++ class) was duplicating memory
  when it was not supposed to. A new implementation avoids copying memory and
  thus makes it faster when trying to fit a larger object.

* Users can define offset terms as function of the vector of sufficient
  statistics.
  
