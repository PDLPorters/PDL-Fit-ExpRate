use strict;
use warnings;

# Totally unnecessary, but needed for to keep Module::Build happy and
# for CPAN indexing. Notice how the VERSION declaration comes before the
# package declaration? Yeah, that's the only way to get it to work. :-/
our $VERSION = "0.03";
package PDL::Fit::ExpRate;

# Invoke PDL::PP
my $base_name;
BEGIN {
	# .PL scripts are sent their filename, sans the .PL part. That's almost what
	# PDL::PP expects to see, so massage it into the proper form:
	$base_name = $ARGV[0];
	$base_name =~ s/\.pm//;
	
	# Handle backslashes for Windows paths:
	$base_name =~ s/\\/\\\\/g;
}
use PDL::PP (__PACKAGE__, __PACKAGE__, $base_name);
# The string passed to pp_setversion must include the quotes to ensure
# that it is processed as a string from the actual module:
pp_setversion("\"$VERSION\"");

# Add the .xs file to the cleanup lists:
use Module::Build;
my $build = Module::Build->current;
$build->add_to_cleanup("$base_name.xs");

pp_addpm <<'DOCUMENTATION';

our $VERSION = "0.03";

=head1 NAME

PDL::Fit::ExpRate - fitting noisy exponential decay

=head1 VERSION

This document describes PDL::Fit::ExpRate version 0.03

=head1 SYNOPSIS

 use PDL;
 use PDL::Fit::ExpRate;
 
 # Load x/y data, or generate some noisy stuff:
 my $xs = sequence(100)/10;
 my ($A, $B, $tau) = (5, 4, -10);
 my $ys = $A + $B * exp($xs / $tau) + $xs->grandom;
 
 # Extract the parameters
 my ($fit_A, $fit_B, $fit_tau) = fit_exp_rate($xs, $ys);
 
 print "A was $A; fit as $fit_A\n";
 print "B was $B; fit as $fit_B\n";
 print "tau was $tau; fit as $fit_tau\n";
 
 # Other useful functions
 
 # Solve for $coefs in $A * $coefs = $y
 # where $A is 3x3
 my $A = pdl q[1 0 0 ; 0 2 0 ; 0 0 3];
 my $y = pdl q[5     ,   -4  ,     6];
 my $coefs = three_by_three_Householder($A, $y);
 
 # Perform a quadratic fit
 # y = coefs[0] + coefs[1] * x + coefs[2] * x**2
 my $coefs = fit_quadratic($xs, $ys);
 
 # Use a callback to get more information from the fit
 ($fit_A, $fit_B, $fit_tau) = fit_exp_rate($xs, $ys,
     run_each_fit => sub {
         my %info = @_;
         print "Full details at end of fit:\n";
         for my $k (sort keys %info) {
             print " $k => $info{$k}\n";
         }
     }
 );

=head1 DESCRIPTION

This module provides a handy function for fitting time series data to a
single exponential growth or decay. A number of methods used by this method
are also made available, including a method for fitting data to a
second-order polynomial, and performing Householder diagonalization on
3x3 matrix systems, and various parameter estimation and single-step functions.

The naive approach to fitting data to an exponential decay involves taking
the logarithm of the y-values and fitting the resulting transformed data to
a line. This has a number of drawbacks. First, it leads to a bias in the
errors: errors in the data associated with smaller y-values receive a
stronger weight than those associated with larger y-values. Second, there is
no simple way to handle an additive offset such as A in C<y = A + B * e**(x/tau)>.
Third and most important, the transformation technique breaks down if the
noise leads to any data values that are negative. The only simple way to
handle such data is to discard it. An approach involving a logarithmic
transformation is biased and can only be used under limited circumstances.

The fitting method supplied by this module, L</fit_exp_rate>, takes a
different and much more robust approach that does not suffer any of those
errors. It performs a least-squares fit to an exponential directly to the
data using an iterative guarded Newton method. The technique gives accurate
results, is highly configurable, and converges rapidly. As an added bonus,
the underlying C functions can be accessed by other PDL modules using
L<PDL::PP::Include::Fit::ExpRate>.

=cut

DOCUMENTATION

# Needed for epsilon
pp_addhdr("#include <float.h>\n");

=head2 FUNCTIONS

=cut

#################################
# 3x3 Household diagonalization #
#################################

pp_addhdr pp_line_numbers(__LINE__, q{
	/* Uses Householder method to diagonalize A*x = y. This modifies A and y
	 * in-place before doing back-substitution to solve for x (which it also
	 * modifies in-place). */ 
	void _three_by_three_Householder (double A[3][3], double y[3], double x[3]) {
		double alpha, gamma;
		int j, k, n, m;
		double v[3];
		double norm_of_max_column, current_norm;
		
		/* This needs to be volatile so that in comparison with zero on line 140
		 * is correct. The better solution would be to perform column pivoting,
		 * but for now this gets me around the -O3 optimization issue. */
		volatile double beta;
		
		/* For column pivotin, we start with the columns in-order. */
		int column_order[3] = {0, 1, 2};
		int real_n;
		
		/* -- Upper-triangulize -- */
		
		for (n = 0; n < 3; n++) {     /* Loop over columns */
			
			/* -- Calculate pivot -- */
			norm_of_max_column = 0;
			for (m = n; m < 3; m++) {
				norm_of_max_column += A[column_order[n]][m]*A[column_order[n]][m];
			}
			for (m = n+1; m < 3; m++) {
				current_norm = 0;
				for (j = n; j < 3; j++) {
					current_norm += A[column_order[m]][j]*A[column_order[m]][j];
				}
				/* pivot if some other column has larger norm */
				if (norm_of_max_column < current_norm) {
					int temp = column_order[n];
					column_order[n] = column_order[m];
					column_order[m] = temp;
				}
			}
			real_n = column_order[n];
			
			/* Compute the Householder vector, v, for the current column; 
			 * alpha is a temporary calculation for this */
			alpha = 0;
			v[0] = v[1] = v[2] = 0;
			for (m = n; m < 3; m++) {
				alpha += A[real_n][m] * A[real_n][m];
				v[m] = A[real_n][m];
			}
			alpha = sqrt(alpha);
			if (A[real_n][real_n] > 0)
				alpha = -alpha;
			v[n] -= alpha;
			
			beta = 0;
			for (m = n; m < 3; m++) {
				beta += v[m] * v[m];
			}
			
			// Only work with values of beta that are definitely nonzero
			if (beta > 0) {

				// apply the Householder transformation to the
				// remaining submatrix, column j, row k
				for (j = n; j < 3; j++) {
					// Get the dot product of v and the jth
					// column of the submatrix, called gamma:
					gamma = 0;
					for (k = n; k < 3; k++) gamma += v[k] * A[column_order[j]][k];
					
					// on gamma, from the submatrix's jth
					// now subtract a scaled version of v, based
					// column:
					for (k = n; k < 3; k++)
						A[column_order[j]][k] -= 2 * gamma / beta * v[k];
				}
				
				// -- apply the Householder transformation to y -- //
				
				// Get the dot product of v and y:
				gamma = 0;
				for (k = n; k < 3; k++) gamma += v[k] * y[k];
				
				// subtract a scaled copy of v from y
				for (k = n; k < 3; k++) {
					y[k] -= (2 * gamma / beta) * v[k];
				}
			}
		}
		
		// -- back-substitution -- //
		
		// here I use alpha as a temporary variable
		
		// Start at the lowest-right diagonal entry:
		for (j = 2; j > -1; j--) {
			alpha = y[j];
			for (k = 2; k > j; k--) {
				// This part of the loop is not entered until x[k]
				// is already defined
				alpha -= A[column_order[k]][j] * x[column_order[k]];
			}
			x[column_order[j]] = alpha / A[column_order[j]][j];
		}
	}
});

=head2 three_by_three_Householder

=cut

my $documentation = <<'DOCUMENTATION';

=for ref

Given A and y, solves A * coef = y

This function performs Householder diagonalization on a 3x3 matrix A with
the given "answer" y. It effectively inverts the matrix, leading to a
solution for coef that can be effectively written as:

 coef = A**-1 * y

where A**-1 is the matrix inverse of A and multiplication is matrix
multiplication.

Here's an exampe of how to use it. I've deliberately made the A matrix
simple so that you can perform the inversion in your head, but A need not
be sparse.

=for example

 my $A = pdl q[1 0 0 ; 0 2 0 ; 0 0 3];
 my $y = pdl q[5     ,   -4  ,     6];
 my $coefs = three_by_three_Householder($A, $y);
 my $expected = pdl(5, -2, 2);
 print "Got $coefs; expected $expected\n";

=cut

DOCUMENTATION

pp_def('three_by_three_Householder',
	Pars => 'A(n=3,m=3); y(n=3); [o] coef(n=3)',
	GenericTypes => ['D'],
	Doc => $documentation,
	Code => pp_line_numbers(__LINE__, q{
		double tmp_A [3][3];
		double tmp_y [3];
		double tmp_c [3];
		threadloop %{
			/* Copy the values into the temporary matrices */
			loop (n) %{
				tmp_y[n] = $y();
				loop (m) %{
					tmp_A[n][m] = $A();
				%}
			%}
			
			/* Call the C-method */
			_three_by_three_Householder(tmp_A, tmp_y, tmp_c);
			loop(n) %{
				$coef() = tmp_c[n];
			%}
		%}
	}),
);

###########################
# Quadratic curve fitting #
###########################
# Uses the 3x3 Householder

pp_addhdr pp_line_numbers(__LINE__, q{
	void _fit_quadratic (double * xs, double * ys, int N_points, double * coefs) {
		double A[3][3];
		double y[3];
		int i, n, m;
		double x, x_sq;
		
		/* Zero-out A and y */
		for (n = 0; n < 3; n++) {
			y[n] = 0;
			for (m = 0; m < 3; m++) {
				A[n][m] = 0;
			}
		}
		
		/* Accumulate A and y */
		for (i = 0; i < N_points; i++) {
			double x_sq = xs[i] * xs[i];
			y[0] += ys[i];
			y[1] += xs[i] * ys[i];
			y[2] += x_sq * ys[i];
			A[0][1] += xs[i];
			A[0][2] += x_sq;
			A[1][2] += x_sq * xs[i];
			A[2][2] += x_sq * x_sq;
		}
		
		/* Fill in the non-unique matrix elements */
		A[0][0] = N_points;
		A[1][0] = A[0][1];
		A[1][1] = A[0][2];
		A[2][0] = A[0][2];
		A[2][1] = A[1][2];
		
		_three_by_three_Householder(A, y, coefs);
	}
});

=head2 fit_quadratic

=cut

$documentation = <<'DOCUMENTATION';

=for ref

Given x and y data, determines A, B, and C for y = A + B x + C x**2

This function determines the least-squares coefficients for a quadratic fit
to a given set of x-y data. This is a simple PDL wrapper for the internal
C function C<exprate_fit_quadratic>. This usage is pretty simple:

=for example

 # Make some noisy data
 my $xs = sequence(50);
 my $ys = 5 + 3 * $x + 4 * $x**2 + grandom($xs);
 
 # should be (5, 3, 4)
 my $coefs = fit_quadratic($xs, $ys);

=cut

DOCUMENTATION

pp_def('fit_quadratic',
	Pars => 'xs(n); ys(n); [o] coefs(m=3)',
	GenericTypes => ['D'],
	Doc => $documentation,
	Code => pp_line_numbers(__LINE__, q{
		_fit_quadratic($P(xs), $P(ys), $SIZE(n), $P(coefs));
	}),
);

################################
# Initial parameter estimation #
################################
# Uses the quadratic fitting

pp_addhdr pp_line_numbers(__LINE__, q{
	double accum_exp_fit_sum_sq_err (double * xs, double * ys, int n_values,
		double A, double B, double lambda
	) {
		int i;
		double total_sq_err = 0, err;
		for (i = 0; i < n_values; i++) {
			err = A + B * exp(lambda * xs[i]) - ys[i];
			total_sq_err += err*err;
		}
		return total_sq_err;
	}
	
	void _exp_fit_estimate_parameters(
		/* Input parameters */
		double * xs, double * ys, int length,
		/* Output containers */
		double * A, double * B, double * lambda, double * sum_sq_error
	) {
		double coefs[3];
		_fit_quadratic(xs, ys, length, coefs);
		
		/* I have found that I tend to get the best results when I use
		 * the average x-position. */
		double x = 0;
		int i = 0;
		for (i = 0; i < length; i++) x += xs[i];
		x /= (double)length;
		
		/* Compute the coefficients by approximating the exponential by
		 * a quadratic fit, evaluated at the average x position. */
		*lambda = 2.0 * coefs[2] / (coefs[1] + 2.0 * coefs[2] * x);
		*B = (coefs[1] + 2.0 * coefs[2] * x)
			/ ((*lambda) * exp(*lambda * x));
		*A = coefs[0] + coefs[1] * x + coefs[2] * x*x
			- (*B) * exp(*lambda * x);
		(*sum_sq_error) = accum_exp_fit_sum_sq_err(xs, ys, length,
			(*A), (*B), (*lambda));
	}
});

=head2 exp_fit_estimate_parameters

=cut

$documentation = <<'DOCUMENTATION';

=for ref

Given x and y data, determines A, B, and lambda for y = A + B exp(x*lambda)

This function obtains an initial estimate for the parameters of an exponential
decay. It also returns the sum of the squared errors for the fit. It obtains
parameter estimates by fitting the data to a quadratic function, then mapping the
coefficients of the quadratic to the coefficients of an exponential function
using a Taylor Series. For example:

=for example

 # Make some noisy data
 my $xs = sequence(50);
 my $ys = 5 + 3 * exp($xs * 0.04) + grandom($xs) * 0.1;
 print "        A      B   lambda    err\n";
 
 # Print the original parameters
 printf "gave: %1.3f  %1.3f  %1.3f\n", 5, 3, 0.04;
 
 # Get estimates for parameters
 my ($A, $B, $lambda, $sum_sq_err) = exp_fit_estimate_parameters($xs, $ys);
 printf "est:  %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;

This is a simple PDL wrapper for the internal C function
C<exprate_estimate_parameters> provided in L<PDL::PP::Include::Fit::ExpRate>.

=cut

DOCUMENTATION

pp_def('exp_fit_estimate_parameters' =>
	Pars => 'xs(n); ys(n); [o] As(); [o] Bs(); [o] lambdas(); [o] sum_sq_errors()',
	GenericTypes => ['D'],
	Code => q{
		/* Really simple: just call the C function */
		_exp_fit_estimate_parameters($P(xs), $P(ys), $SIZE(n),
				$P(As), $P(Bs), $P(lambdas), $P(sum_sq_errors));
	},
	Doc => $documentation,
);

#############################
# Single Newton-method step #
#############################
# Uses 3x3 Householder

pp_addhdr pp_line_numbers(__LINE__, q{
	int _exp_fit_newton_method_step (
		/* Input variables */
		double * xs, double * ys, int length, double trust_radius,
		/* Output variables */
		double * A, double * B, double * lambda, double * sum_sq_error
	) {
		double Hessian[3][3], neg_gradient[3], coefs[3], step[3];
		int i;
		
		/**********************/
		/* Newton-method step */
		/**********************/
		
		/* Compute the gradient and the Hessian */
		Hessian[0][1] = Hessian[0][2] = Hessian[1][1]
			= Hessian[1][2] = Hessian[2][2]
			= neg_gradient[0] = neg_gradient[1] = neg_gradient[2] = 0;
		for (i = 0; i < length; i++) {
			double x = xs[i];
			double the_exp = exp(x * (*lambda));
			double dy = (*A) + (*B) * the_exp - ys[i];
			
			/* Terms for the negative neg_gradient */
			neg_gradient[0] -= dy;
			neg_gradient[1] -= dy * the_exp;
			neg_gradient[2] -= dy * x * (*B) * the_exp;
			
			/* Terms for the Hessian */
			double weird = (*B) * the_exp + dy;
			Hessian[0][1] += the_exp;
			Hessian[0][2] += x * (*B) * the_exp;
			Hessian[1][1] += the_exp * the_exp;
			Hessian[1][2] += x * the_exp * weird;
			Hessian[2][2] += x * x * (*B) * the_exp * weird;
		}
		
		/* finish the (symmetric) Hessian */
		Hessian[0][0] = length;
		Hessian[1][0] = Hessian[0][1];
		Hessian[2][0] = Hessian[0][2];
		Hessian[2][1] = Hessian[1][2];
		
		/* Compute the step size using Householder diagonalization */
		_three_by_three_Householder(Hessian, neg_gradient, step);
		
		/**********************/
		/* Gaurded correction */
		/**********************/
		
		/* Gaurd the step so it's not bigger than the trust radius */
		double correction = 1.0;
		double potential_correction;
		/* check A step */
		potential_correction = trust_radius * fabs(*A) / fabs(step[0]);
		if (potential_correction < correction) correction = potential_correction;
		/* check B step */
		potential_correction = trust_radius * fabs(*B) / fabs(step[1]);
		if (potential_correction < correction) correction = potential_correction;
		/* check lambda step */
		potential_correction = trust_radius * fabs(*lambda) / fabs(step[2]);
		if (potential_correction < correction) correction = potential_correction;
		
		/* Take the recommended step, with applicable correction */
		(*A) += step[0] * correction;
		(*B) += step[1] * correction;
		(*lambda) += step[2] * correction;
		
		/* Compute the new score */
		*sum_sq_error = accum_exp_fit_sum_sq_err(xs, ys, length, *A, *B, *lambda);
		
		return (correction < 1.0);
	}
});

=head2 exp_fit_newton_method_step

=cut

$documentation = <<'DOCUMENTATION';

=for ref

Takes a single gaurded Newton-method step for the iterative fitting process

While the above method, C<exp_fit_estimate_parameters>, can give reasonable
estimates for an exponential fit, a proper fit involves an iterative Newton
method. Such a method evaluates the local gradient of the error and takes a
step in the direction toward the estimated minimum of the error. The method
can only estimate where the minimum is located, so it limits the size of its
steps so as not to go too far. While the full fittinig routine performs many
steps, this function performs only a single step.

This function takes a number of parameters, and returns quite a few as well.
The input values include the data as well as the previous best guess at the
parameters for A, B, and lambda in the equation

 y = A + B * exp(x * lambda)

The function returns the result of a single gaurded Newton method iteration.
Such a step is guaranteed to fall within the specified trust radius, which I
explain more in a moment. The function also returns the sum of the squared
errors for the new parameters.

=for example

 # Make some noisy data
 my $xs = sequence(50);
 my $ys = 5 + 3 * exp($xs * 0.04) + grandom($xs) * 0.1;
 print "          A      B   lambda  err\n";
 
 # Get initial estimates for parameters
 my ($A, $B, $lambda, $sum_sq_err) = exp_fit_estimate_parameters($xs, $ys);
 printf "step 0: %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 
 # Perform many single steps to see how it goes
 my $trust_radius = 0.1;  # trust radius of 10%
 for (1 .. 9) {
   ($A, $B, $lambda, $sum_sq_err)
     = exp_fit_newton_method_step($xs, $ys, $trust_radius, $A, $B, $lambda);
   printf "step $_: %1.3f  %1.3f  %1.3f  %1.3f\n", $A, $B, $lambda, $sum_sq_err;
 }
 
 # Print the original parameters
 print '=' x 35, "\n";
 printf "gave:   %1.3f  %1.3f  %1.3f\n", 5, 3, 0.04;

The trust radius ensures that each of the parameters does not change by more
than the specified amount in a given step. The basic idea is that we want to
find the minimum in the surface representing the sum of the squared errors,
and we do this by locally approximating the surface by a quadratic. This local
approximation is only accurate within a small parameter radius, so if the
minimum of the local quadratic is too far away, it would be unwise to go all
the way to that location. A trust radius of 0.1 (i.e. 10%) means you do not
want your estimates to change by more than 10%.

You should avoid a trust radius larger than 1. Such a large trust radius would
make it possible to change the sign of your parameters. If you find that you
are getting the wrong sign for your parameters, you should find a better initial
estimate rather than increase your trust radius.

=cut

DOCUMENTATION


pp_def('exp_fit_newton_method_step' =>
	Pars => 'xs(n); ys(n); trust_radii(); in_As(); in_Bs(); in_lambdas();
			[o] out_As(); [o] out_Bs(); [o] out_lambdas(); [o] sum_sq_err()',
	GenericTypes => ['D'],
	Doc => $documentation,
	Code => q{
		/* The A/B/lambda values are modified in-place by the original function
		 * so copy the input values to the output values before calling the
		 * function. */
		$out_As() = $in_As();
		$out_Bs() = $in_Bs();
		$out_lambdas() = $in_lambdas();
		
		_exp_fit_newton_method_step($P(xs), $P(ys), $SIZE(n), $trust_radii(),
			&($out_As()), &($out_Bs()), &($out_lambdas()), &($sum_sq_err()));
	},
);

################################
# Code to handle the callbacks #
################################
# I've broken these out into #defines because I wanted to clarify the
# actual fit_exp_rate code, and to avoid code duplication on
# PER_ITERATION_CALLBACK
pp_addhdr pp_line_numbers(__LINE__, q{
	#define DO_ITERATION_CALLBACK do { \
		dSP; \
		ENTER; \
		SAVETMPS; \
		\
		PUSHMARK(SP); \
		EXTEND(SP, 14); \
		PUSHs(sv_2mortal(newSVpv("A", 1))); \
		PUSHs(sv_2mortal(newSVnv(A))); \
		PUSHs(sv_2mortal(newSVpv("B", 1))); \
		PUSHs(sv_2mortal(newSVnv(B))); \
		PUSHs(sv_2mortal(newSVpv("lambda", 6))); \
		PUSHs(sv_2mortal(newSVnv(lambda))); \
		PUSHs(sv_2mortal(newSVpv("round", 5))); \
		PUSHs(sv_2mortal(newSViv(counter))); \
		PUSHs(sv_2mortal(newSVpv("sum_sq_err", 10))); \
		PUSHs(sv_2mortal(newSVnv(sum_sq_error))); \
		PUSHs(sv_2mortal(newSVpv("old_sum_sq_err", 14))); \
		PUSHs(sv_2mortal(newSVnv(old_sum_sq_error))); \
		PUSHs(sv_2mortal(newSVpv("threshold", 9))); \
		PUSHs(sv_2mortal(newSVnv(threshold))); \
		if (force_next_round) { \
			XPUSHs(sv_2mortal(newSVpv("forced_round", 12))); \
			XPUSHs(sv_2mortal(newSViv(force_next_round))); \
		} \
		PUTBACK; \
		\
		call_sv(*per_iteration_code_ref_p, G_DISCARD); \
		\
		FREETMPS; \
		LEAVE; \
	} while (0)
	
	#define DO_FIT_CALLBACK do { \
		dSP; \
		ENTER; \
		SAVETMPS; \
		\
		PUSHMARK(SP); \
		EXTEND(SP, 20); \
		PUSHs(sv_2mortal(newSVpv("fit_count", 9))); \
		PUSHs(sv_2mortal(newSViv(total_count))); \
		PUSHs(sv_2mortal(newSVpv("N_fits", 6))); \
		PUSHs(sv_2mortal(newSViv(total_expected))); \
		PUSHs(sv_2mortal(newSVpv("sum_sq_err", 10))); \
		PUSHs(sv_2mortal(newSVnv(sum_sq_error))); \
		PUSHs(sv_2mortal(newSVpv("old_sum_sq_err", 14))); \
		PUSHs(sv_2mortal(newSVnv(old_sum_sq_error))); \
		PUSHs(sv_2mortal(newSVpv("threshold", 9))); \
		PUSHs(sv_2mortal(newSVnv(threshold))); \
		PUSHs(sv_2mortal(newSVpv("N_rounds", 8))); \
		PUSHs(sv_2mortal(newSViv(counter))); \
		PUSHs(sv_2mortal(newSVpv("A", 1))); \
		PUSHs(sv_2mortal(newSVnv(A))); \
		PUSHs(sv_2mortal(newSVpv("B", 1))); \
		PUSHs(sv_2mortal(newSVnv(B))); \
		PUSHs(sv_2mortal(newSVpv("tau", 3))); \
		PUSHs(sv_2mortal(newSVnv(1.0/lambda))); \
		PUSHs(sv_2mortal(newSVpv("converged", 9))); \
		PUSHs(sv_2mortal(newSVnv(converged))); \
		PUTBACK; \
		\
		callback_return = call_sv(*per_fit_code_ref_p, G_SCALAR); \
		\
		SPAGAIN; \
		\
		if (callback_return != 1) croak ("How did I get more than one arg returned?"); \
		\
		/* See if they want to quit early */ \
		callback_return = POPi; \
		if (callback_return == 0) goto QUIT; \
		\
		PUTBACK; \
		FREETMPS; \
		LEAVE; \
	} while (0)
});


#############################################
# The main all-encompasing fitting function #
#############################################

=head2 fit_exp_rate

=cut

$documentation = <<'DOCUMENTATION';

=for ref

Fits your data to an exponential decay or growth

This is the primary method of PDL::Fit::ExpRate. This method takes your
x and y data and computes a best-fit for the parameters A, B, and tau in
the equation

                  x[i] / tau
 y[i] = A + B * e

If a fit fails to converge, the values of A, B, and tau for that dataset
will be set to C<BAD> (see L<PDL::Bad> for details).

In additon to the x/y data, you can optionally include A, B, and tau.
For most PDL methods, this provides a means to pre-allocate (or use
previously allocated) memory for the output, which can be important for
large data operations. That's not likely to be the case for this method,
but the mechanism is provided for consistency. This can also be used if
you want to provide your own initial estimates for A, B, and tau, which
you may find useful if the quadratic estimate method is not working well
for your data. The method assumes the latter, but if you are simply want
to pass along pre-allocated piddles, see the C<pre_allocated> flag below.
Either way, the method will overwrite the data in A, B, and tau, as they
are considered output piddles.

After piddle arguments x and y (and optionally A, B, and tau), you can
specify key/value pairs that govern how the fit is performed. Options
include:

=over

=item pre_allocated

If you pre-allocated piddles for A, B, and tau, and are not providing
them as initial estimates, then you should include this key with a
boolean true value, i.e. 1.

=item trust_radius

The trust radius for the fit(s), which governs how much the parameters
will evolve with a single iteration. See the documentation for
C<exp_fit_newton_method_step> for a more complete explanation of this
parameter. The default value is 0.1, i.e. 10%.

=item iterations

The maximum number of iterations, i.e. the fit must converge before this
number of iterations or it will be considered a failed fit. The default
value is 50.

=item threshold

The sum of the squared errors will vary as the fit iterates. As the fit
converges, the deviations will get very small. The C<threshold> indicates
the deviations that are acceptable to declare convergence. This is not the
percent difference in the parameters, but in the sum of the squared errors.
The default is 0.001.

Note that for noiseless data, the fitting routine may run into issues with
bit noise in its calculations for the sum of the squared errors. If it
detects that the sum of the squared errors has fallen into bitnoise, it
will deem the system converged and disregard the threshold.

=item min_lambda, max_lambda

You may have physical reasons for expecting your exponential rates to fall
within certain ranges. If so, you can provide the min and max values of
lambda that are allowed. (Lambda is the inverse of tau.) The magnitudes
of these numbers will be compared with the magnitudes of the lambda value
with each iteration of the fit, and if the fit value falls outside of the
given range, the fit will be considered a failure. If you do not want to
impose any limitations, set both of these to zero. The default value for
max_lambda is zero (no limit) but the default value for min_lambda is
1e-8.

=item run_each_fit

A callback that gets called at the completion of each fit. This was
implemented primarily for use in progress indicators and GUI applications.
However, it is also useful for getting more information about each fit
than is returned by the method itself. The items passed to the callback are
key/value pairs (i.e. something that you can assign to a hash) that
include:

  N_fits     - number of fits we're doing
  fit_count  - number of fits we've completed
  A          - final value for A for this fit
  B          - final value for B for this fit
  tau        - final value for tau for this fit
  converged  - boolean indicating if the routine will declare
               that this fit has converged
  sum_sq_err - last sum of the squared deviations for the
               just-completed fit
  old_sum_sq_err
             - previous sum of the squared deviations for the
               just-completed fit
  threshold  - threshold, see option above
  N_rounds   - number of rounds required for convergence

Note that if C<converged> is false, the method will eventually return
values for A, B, and tau that are BAD. This callback provides a way to
get the last (but unsuccessful) values of A, B, and tau, which may be
helpful in certain circumstances. Usually, converged will be true if
(sum_sq_err - old_sum_sq_err) / old_sum_sq_err < threshold. It might also
be true if the data is noise-free and sum_sq_err has fallen to the level
of bit noise.

The return value of this callback is meaningful. If it is false, the
PDL function will abort without threading over any more datasets. If it
is true, the PDL function will continue. This has seen use in GUI
programs that initiate a large number of fits. The callback was used to
keep the GUI responsive, and to make it possible for the user to cancel
the fits mid-calculation.

A simple progress indicator might look like this:

 my ($As, $Bs, $taus) = fit_exp_rate($xs, $ys,
     run_each_fit => sub {
         my %info = @_;
         print "Completed $info{fit_count} of $info{N_fits}\n";
         return 1;
     },
 );

=item run_each_iteration

A callback that gets called at the completion of each iteration of each
fit. As this will get called with each round of the iteration, an
iteration callback can have a significant performance impact. As such,
this is primarily useful as a debugging aid when your fit fails to
converge. In such circumstances, you can add an iteration callback to
print how the parameters are evolving, which may help figure out how to
tweak the parameters that govern the iterative process.

The items passed to the callback are key/value pairs (i.e. something
that you can assign to a hash) that include:

  A           - current value for A for this fit
  B           - current value for B for this fit
  lambda      - current value for 1/tau for this fit
  round       - iteration count
  sum_sq_err  - last sum of the squared deviations for the
                just-completed fit
  old_sum_sq_err
              - previous sum of the squared deviations for
                the just-completed fit
  threshold   - threshold, see option above
  force_round - boolean indicating whether the next round
                is forced, which happens if the Newton
                method step was confined by the trust radius

Of these, C<force_next_round> probably needs a bit of clarification. If
the Newton method step was truncated by the C<trust_radius>, then the
change in the sum of the squared errors may not have changed much. In
that case, the threshold limit could erroneously indicate convergence.
For this reason, the method will always perform another iteration when
the step gets truncated by the trust radius. The value of C<force_next_round>
will indicate if this is the case.

=back

All of these values are either scalars or function references, but are
not piddles. So, unlike the trust radius in C<exp_fit_newton_method_step>,
if you utilize PDL's threading to perform exponential fits over multiple
datasets in one call, they will all use this same value of the trust_radius
and the other parameters.

Notice that, unlike exp_fit_newton_method_step, this method does not
return the sum of the squared errors. If you want that information, you
should supply a per-fit callback which collects the result.

Here are some examples of how to use the function:

 # Make some noisy data
 my $xs = sequence(500);
 my $ys = 5 + 3 * exp($xs * -0.05) + grandom($xs) * 0.5;
 print "   A      B      tau\n";
 
 # Use default parameters
 my ($A, $B, $tau) = fit_exp_rate($xs, $ys);
 printf " %1.3f  %1.3f  %1.3f  -> defaults\n", $A->sclr, $B->sclr, $tau->sclr;
 
 # See how things as we tighten the threshold a bit
 for my $threshold (0.1, 0.01, 0.001, 0.0001, 0.00001) {
     ($A, $B, $tau) = fit_exp_rate($xs, $ys,
       threshold => $threshold
     );
     printf " %1.3f  %1.3f  %1.3f  -> threshold = $threshold\n", $A->sclr, $B->sclr, $tau->sclr;
 }
 
 print '=' x 35, "\n";
 printf " %1.3f  %1.3f  %1.3f  -> values for generating data\n", 5, 3, -1/0.05;
 
 # Use callback to track the convergence of A
 my @As;
 ($A, $B, $tau) = fit_exp_rate($xs, $ys,
     run_each_iteration => sub {
         my %info = @_;
         push @As, $info{A};
     }
 );
 print "\nA:   ", join ("\n  -> ", @As), "\n";
 
 # Use callback to print diagnostics
 print '-' x 30, "\n";
 ($A, $B, $tau) = fit_exp_rate($xs, $ys,
    iterations => 100,
	run_each_iteration => sub {
		my %info = @_;
		return if $info{round} % 5 != 0;
		print "Round $info{round}:\n";
		delete $info{round};
		for my $k (sort keys %info) {
			print "  $k: $info{$k}\n";
		}
	},
 );

=cut

DOCUMENTATION

pp_def('fit_exp_rate',
	Pars => 'xs(n); ys(n); [o] As(); [o] Bs(); [o] taus(); int [o] is_bad()',
	GenericTypes => ['D'],
	OtherPars => 'SV * options_sv',
	PMCode => pp_line_numbers(__LINE__, q{
		use Carp 'croak';
		use PDL;
		use strict;
		use warnings;
		
		sub PDL::fit_exp_rate {
			# Load the data:
			my ($xs, $ys) = (shift, shift);
			$xs = PDL::Core::topdl($xs);
			$ys = PDL::Core::topdl($ys);
			
			my ($A, $B, $tau, $provided_outputs);
			if (eval {$_[0]->isa('PDL')} ) {
				($A, $B, $tau) = (shift, shift, shift);
				$provided_outputs = 1;
			}
			else {
				my @dims = $xs->dims;
				shift @dims;
				my @y_dims = $ys->dims;
				shift @y_dims;
				for (my $i = 0; $i < @y_dims; $i++) {
					$dims[$i] = $y_dims[$i]
						if not defined $dims[$i]
							or $dims[$i] == 1 and $y_dims[$i] > 1;
				}
				
				# Handle single-input case
				@dims = (1) unless @dims;
				
				$A = zeroes(@dims);
				$B = zeroes(@dims);
				$tau = zeroes(@dims);
			}
			my $is_bad = $A->ones;
			
			croak ("You must provide key => value pairs")
				unless @_ % 2 == 0;
			
			my %opts = @_;
			
			# If A, B, and tau were not provided, then set the
			# pre-allocated flag.
			$opts{pre_allocated} = 1 if not $provided_outputs;
			# Make it integer 0 or 1 to reflect boolean state
			$opts{pre_allocated} = $opts{pre_allocated} ? 1 : 0;
			
			# Determine the total number of rounds for the run
			my $n_rounds = 1;
			$n_rounds *= $_ for ($A->dims);
			$opts{N_rounds} = $n_rounds;
			
			PDL::_fit_exp_rate_int($xs, $ys, $A, $B, $tau, $is_bad, \%opts);
			
			# If they exited early, not all the return values will be good
			if (any $is_bad) {
				$A = $A->setbadif($is_bad);
				$B = $B->setbadif($is_bad);
				$tau = $tau->setbadif($is_bad);
			}
			return ($A, $B, $tau);
		}
	}),
	Code => pp_line_numbers(__LINE__, q{
		
		double step[3];
		double A, B, lambda, old_sum_sq_error, sum_sq_error;
		
		/* If we take a guarded step, then we force the next round. This is so
		 * that the fit won't erroneously think it's converged just because the
		 * guarded step led to a very small change in the sum_sq_error */
		int force_next_round, converged;
		
		/**************************************/
		/* Process all the optional arguments */
		/**************************************/
		
		HV * options_hv = (HV*)SvRV($COMP(options_sv));
		
		/* Get the per-fit and per-round callbacks, if given */
		SV ** per_fit_code_ref_p = hv_fetchs(options_hv, "run_each_fit", 0);
		SV ** per_iteration_code_ref_p = hv_fetchs(options_hv, "run_each_iteration", 0);
		int callback_return;
		
		/* Minimum lambda before failing */
		SV ** min_lambda_p = hv_fetchs(options_hv, "min_lambda", 0);
		double min_lambda = min_lambda_p ? SvNV(*min_lambda_p) : 1e-8;
		
		/* Maximum lambda before failing; value 0 means any value is allowed */
		SV ** max_lambda_p = hv_fetchs(options_hv, "max_lambda", 0);
		double max_lambda = max_lambda_p ? SvNV(*max_lambda_p) : 0;
		
		/* Get the trust radius or choose a default of 10% */
		SV ** trust_radius_p = hv_fetchs(options_hv, "trust_radius", 0);
		double trust_radius = trust_radius_p ? SvNV(*trust_radius_p) : 0.1;
		
		/* Set the total number of expected rounds */
		SV ** n_rounds_p = hv_fetchs(options_hv, "N_rounds", 0);
		int total_expected = n_rounds_p ? SvIV(*n_rounds_p) : 1;
		int total_count = 0;
		
		/* Pull max_iterations from the hash or set to a reasonable default */
		SV ** max_it_p = hv_fetchs(options_hv, "iterations", 0);
		int max_iterations = max_it_p ? SvIV(*max_it_p) : 50;
		int counter;
		
		/* Pull the threshold from the hash or set to a reasonable default */
		SV ** threshold_p = hv_fetchs(options_hv, "threshold", 0);
		double threshold = threshold_p ? SvNV(*threshold_p) : 0.001;
		
		/* Determine if they provided initial estimates */
		SV ** pre_allocated_p = hv_fetchs(options_hv, "pre_allocated", 0);
		int already_estimated = 1 - SvIV(*pre_allocated_p);
		
		/* If the sum of the squared error drops below this, we're at the
		 * level of bitnoise, so we'll declare victory (convergence). */
		double minimum_calculable_sum_sq_error
			= $SIZE(n)*$SIZE(n) * DBL_EPSILON * DBL_EPSILON;
		/* This value was determined semi-empirically. I feel like the
		 * factor of n**2 should just be n, and I should have a fudge
		 * factor, but the fudge factor seemed to scale with the system
		 * size so it just became n**2. */
		
		/**************/
		/* Threadloop */
		/**************/
		
		threadloop %{
			total_count++;
			
			/* Perform the initial parameter estimate */
			if (already_estimated) {
				A = $As();
				B = $Bs();
				lambda = 1 / $taus();
			}
			else {
				_exp_fit_estimate_parameters($P(xs), $P(ys), $SIZE(n), 
					&A, &B, &lambda, &sum_sq_error);
			}
			
			/* Make sure we enter the loop */
			force_next_round = 1;
			counter = 0;
			
			/**************************************************/
			/* Iterative Newton method to best-fit parameters */
			/**************************************************/
			
			while(
				force_next_round || 
				sum_sq_error > minimum_calculable_sum_sq_error &&
				fabs(old_sum_sq_error - sum_sq_error) > old_sum_sq_error * threshold
			) {
				counter++;
				if (counter > max_iterations) break;
				if (fabs(lambda) < min_lambda) break;
				if (max_lambda > 0 && fabs(lambda) > max_lambda) break;
				
				/* Per-iteration callback */
				if (per_iteration_code_ref_p != NULL) DO_ITERATION_CALLBACK;
				
				/* Call the single stepper */
				old_sum_sq_error = sum_sq_error;
				force_next_round = _exp_fit_newton_method_step(
					$P(xs), $P(ys), $SIZE(n), trust_radius,
					&A, &B, &lambda, &sum_sq_error);
			}

			/* Final per-iteration callback */
			if (per_iteration_code_ref_p != NULL) DO_ITERATION_CALLBACK;

			/********************************/
			/* Store the final coefficients */
			/********************************/
			
			$As() = A;
			$Bs() = B;
			$taus() = 1/lambda;
			$is_bad() = fabs(lambda) < min_lambda
					|| (max_lambda > 0 && fabs(lambda) > max_lambda)
					|| counter > max_iterations;
			converged = !$is_bad();
			
			/* Per-fit callback */
			if (per_fit_code_ref_p != NULL) DO_FIT_CALLBACK;
		%}
	
	QUIT: ;
	
	}),
);

#####################################
# Low-level function pointer copyer #
#####################################

pp_add_boot pp_line_numbers(__LINE__, q{
	
	/* Set global variables with the function addresses, used (potentially)
	 * by PDL::PP::Include::Fit::ExpRate-based code */
	
	/* householder */
	SV * address_sv = get_sv("PDL::Fit::ExpRate::__householder_func_addr", GV_ADD);
	sv_setiv(address_sv, PTR2IV(&_three_by_three_Householder));
	
	/* quadratic fit */
	address_sv = get_sv("PDL::Fit::ExpRate::__quadratic_func_addr", GV_ADD);
	sv_setiv(address_sv, PTR2IV(&_fit_quadratic));
	
	/* exponential error accumulator */
	address_sv = get_sv("PDL::Fit::ExpRate::__accum_func_addr", GV_ADD);
	sv_setiv(address_sv, PTR2IV(&accum_exp_fit_sum_sq_err));
	
	/* parameter estimator */
	address_sv = get_sv("PDL::Fit::ExpRate::__estimate_func_addr", GV_ADD);
	sv_setiv(address_sv, PTR2IV(&_exp_fit_estimate_parameters));
	
	/* newton step */
	address_sv = get_sv("PDL::Fit::ExpRate::__newton_func_addr", GV_ADD);
	sv_setiv(address_sv, PTR2IV(&_exp_fit_newton_method_step));
});

pp_addpm <<'DOCUMENTATION';

=head1 SEE ALSO

There are a number of PDL modules for fitting data, including
L<PDL::Fit::Linfit> for linear fitting and L<PDL::Fit::Polynomial> for
fitting data to polynomials. Both of those modules depend on Slatec (though
they may be installed on your system even if you don't have Slatec available).
L<PDL::Fit::LM> provides Levenberg-Marquardt fitting, and
L<PDL::Fit::Gaussian> provides methods for fitting a set of data to a
normal (Gaussian) distribution.

=head1 AUTHOR

David Mertens <dcmertens.perl@gmail.com>

=head1 LICENCE AND COPYRIGHT

Documentation is copyright (c) 2012, David Mertens, all rights reserved.

Code is copyright (c) 2012, Northwestern University, all rights reserved.

This module is free software; you can redistribute it and/or modify it under
the same terms as Perl itself. See L<perlartistic>.

=cut

DOCUMENTATION

pp_done;
