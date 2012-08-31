package PDL::PP::Include::Fit::ExpRate;
use strict;
use warnings;

=head1 NAME

PDL::PP::Include::Fit::ExpRate - build-time includes to expose ExpRate C functions

=head1 SYNOPSIS

 # In your PDL::PP file:
 use PDL::PP::Include::Fit::ExpRate;
 
 # Later in your PP code
 pp_def ('your_pp_function', Pars => 'xs(n); ys(n); ...',
     Code => q{
         ...
         double A, B, lambda, sum_sq_error;
         
         /* Get the initial parameter estimate for an exponential fit */
         exprate_estimate_parameters($P(xs), $P(ys), $SIZE(n), 
             &A, &B, &lambda, &sum_sq_error);
         ...
     },

=head1 DESCRIPTION

L<PDL::Fit::ExpRate> provides a method for fitting a noisy time series
to an exponential decay. The method is built upon a number of simpler
C functions which themselves have PDL methods. This module is meant to
be included in your L<PDL::PP>-based files (i.e. .pd files) so that your
own C code can access the same C functions used by L<PDL::Fit::ExpRate>.

This module accomplishes all of this by calling a number of pp-functions
when you say C<use PDL::PP::Include::Fit::ExpRate>. These pp-functions
add code to the generated XS and Perl-module files generated by PDL::PP.
Having used this module, you can safely use the following C functions
in the code for your C<pp_addxs> and C<pp_def> declarations:

=over

=item exprate_three_by_three_Householder

=for Sig

  Sig: void exprate_three_by_three_Householder (double A[3][3], double y[3], double x[3]);

=item exprate_quadratic_fit

    void exprate_quadratic_fit (double *, double *, int, double *);

=item exprate_accum_sum_sq_err

    double exprate_accum_sum_sq_err (double *, double *, int, double, double, double);

=item exprate_estimate_parameters

    void exprate_estimate_parameters (double *, double *, int, double *, double *, double *, double *);

=item exprate_newton_method_step

    int exprate_newton_method_step (double *, double *, int, double, double *, double *, double *, double *);

=back

=cut

# Make sure they already included PDL::PP
croak("You must include PDL::PP before PDL::PP::Include::Fit::Exprate")
	unless defined $PDL::PP::VERSION;

sub import {

#############################
# Add the function pointers #
#############################
pp_addhdr q{
    void (*exprate_three_by_three_Householder) (double A[3][3], double y[3], double x[3]);
    void (*exprate_quadratic_fit) (double *, double *, int, double *);
    double (*exprate_accum_sum_sq_err) (double *, double *, int, double, double, double);
    void (*exprate_estimate_parameters) (double *, double *, int, double *, double *, double *, double *);
    int (*exprate_newton_method_step) (double *, double *, int, double, double *, double *, double *, double *);
};


##############################################
# XS code that assigns the function pointers #
##############################################
pp_addxs <<XSCODE;

void
__set_up_exprate_pointers()
	CODE:
		/* Assign the pointers */
        exprate_three_by_three_Householder
            = INT2PTR(SvIV(get_sv("PDL::Fit::ExpRate::__householder_func_addr", 0)));
        exprate_quadratic_fit
            = INT2PTR(SvIV(get_sv("PDL::Fit::ExpRate::__quadratic_func_addr", 0)));
        exprate_accum_sum_sq_err
            = INT2PTR(SvIV(get_sv("PDL::Fit::ExpRate::__accum_func_addr", 0)));
        exprate_estimate_parameters
            = INT2PTR(SvIV(get_sv("PDL::Fit::ExpRate::__estimate_func_addr", 0)));
        exprate_newton_method_step
            = INT2PTR(SvIV(get_sv("PDL::Fit::ExpRate::__newton_func_addr", 0)));

XSCODE


############################################
# Module code that performs the assignment #
############################################
pp_addpm q{
    # Include the ExpRate module
    use PDL::Fit::ExpRate;
    
    # Add code to the module that executes the pointer assignment
    __set_up_exprate_pointers();
};

}

1;
