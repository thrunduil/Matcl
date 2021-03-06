#include <stdio.h>

/*! @file input_error_dist_dist.c
 * \brief Error handler for input parameters.
 *
 * <pre>
 * -- SuperLU routine (version 4.4) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * November 20, 2012
 * </pre>
 */

/*! \brief
 *
 * <pre>
 * Purpose   
 * =======   
 *
 * INPUT_ERROR is called if an input parameter has an   
 * invalid value.  A message is printed and execution stops.   
 *
 * Arguments   
 * =========   
 *
 * srname  (input) character*6
 *         The name of the routine which called INPUT_ERROR.
 *
 * info    (input) int
 *         The position of the invalid parameter in the parameter list   
 *         of the calling routine.
 *
 * </pre>
 */
int input_error_dist(char *srname, int *info)
{
    printf("** On entry to %6s, parameter number %2d had an illegal value\n",
		srname, *info);
    return 0;
}
