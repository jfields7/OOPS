#ifndef TYPES_H
#define TYPES_H

// A general enumerator to describe particular outcomes.
enum Result{
  SUCCESS,
  FAILURE,
  BAD_ALLOC,
  OUT_OF_BOUNDS,
  INVALID_STAGE,
  UNRECOGNIZED_PARAMS,
};

enum Boundary{
  NONE,
  LEFT,
  RIGHT,
};

/*const char *ERROR_CODES[] = {
  "Success",
  "Failure",
  "Bad allocation",
  "Out of bounds"
};*/

#endif
