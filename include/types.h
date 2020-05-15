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
  FIELD_EXISTS,
  UNRECOGNIZED_FIELD,
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

// If we're using g++, we have access to its restrict keyword, which can be useful
// for optimization. To keep the code from breaking, we still define it for other
// compilers even though it won't do anything.
#ifndef __GNUG__
#define RESTRICT
#else
#define RESTRICT __restrict__
#endif

#endif
